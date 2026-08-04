## Multiply imputed designs.
##
## The strategy is to make the pooled result *look like* a fitted model: everything downstream reads
## a model through coef(), vcov(), terms() and family(), so a shim answering those four generics lets
## the whole single-design machinery -- contrasts, modifiers, tests, methods -- work unchanged. The
## one genuine difference is the degrees of freedom, which vary from curve point to curve point and
## are threaded through separately.

is_mi_design <- function(design) {
  inherits(design, "svyimputationList") || inherits(design, "DBimputationList")
}

mi_components <- function(design) {
  comps <- design$designs
  if (!is.list(comps) || !length(comps)) {
    stop_svyrcs("the imputed design contains no component designs")
  }
  if (length(comps) < 2L) {
    stop_svyrcs("multiple imputation needs at least 2 imputations; the design has ", length(comps),
                ". With a single completed dataset, pass that design directly.")
  }
  comps
}

## Every imputation is the same sample, so the complete-data degrees of freedom must agree.
shared_degf <- function(comps) {
  dfs <- vapply(comps, survey::degf, numeric(1L))
  if (length(unique(round(dfs, 8))) > 1L) {
    stop_svyrcs("the imputed designs have different degrees of freedom (",
                paste(unique(dfs), collapse = ", "), "). Every imputation must describe the same ",
                "sample, so this usually means the component designs were built differently or ",
                "subset differently.")
  }
  dfs[[1L]]
}

#' Combine fits across imputations by Rubin's rules
#'
#' Internal. Returns an object that behaves enough like a fitted model for the rest of the package to
#' use it, while also carrying the within- and between-imputation covariance matrices that the
#' degrees of freedom need.
#'
#' @param fits List of fitted models, one per imputation.
#' @param degf Complete-data design degrees of freedom.
#'
#' @return An object of class `svyrcs_mifit`, inheriting the class of the component fits.
#'
#' @keywords internal
#' @noRd
pool_fits <- function(fits, degf) {
  m <- length(fits)
  if (m < 2L) stop_svyrcs("pooling needs at least 2 imputations; got ", m)

  ## Compare the coefficient names before building the matrix. vapply() with a fixed FUN.VALUE
  ## length fails first otherwise, with an internal message about lengths that says nothing about
  ## the actual cause.
  coefs <- lapply(fits, stats::coef)
  nms <- names(coefs[[1L]])
  differs <- vapply(coefs, function(b) !identical(names(b), nms), logical(1L))
  if (any(differs)) {
    culprit <- which(differs)[1L]
    stop_svyrcs("the imputations produced models with different coefficients: imputation ", culprit,
                " has ", length(coefs[[culprit]]), " where imputation 1 has ", length(nms),
                ". This happens when a factor has levels present in some imputations but not ",
                "others, or when a fit is rank deficient in only some of them.")
  }

  betas <- matrix(unlist(coefs, use.names = FALSE), nrow = length(nms), dimnames = list(nms, NULL))

  beta_bar <- rowMeans(betas)
  Ubar <- Reduce(`+`, lapply(fits, stats::vcov)) / m
  ## cov() uses the (m - 1) denominator, matching mitools::MIcombine().
  B <- stats::cov(t(betas))
  dimnames(B) <- dimnames(Ubar)
  Tm <- Ubar + (1 + 1 / m) * B

  structure(
    list(coefficients = beta_bar, variance = Tm, Ubar = Ubar, B = B, m = m,
         degf = degf, fits = fits),
    ## Inheriting the component class is what makes effect_measure() and family() dispatch work;
    ## the methods below take precedence over the inherited ones.
    class = c("svyrcs_mifit", class(fits[[1L]]))
  )
}

#' @export
coef.svyrcs_mifit <- function(object, ...) object$coefficients

#' @export
vcov.svyrcs_mifit <- function(object, ...) object$variance

#' @export
terms.svyrcs_mifit <- function(x, ...) stats::terms(x$fits[[1L]], ...)

#' @export
family.svyrcs_mifit <- function(object, ...) stats::family(object$fits[[1L]], ...)

#' @export
model.frame.svyrcs_mifit <- function(formula, ...) stats::model.frame(formula$fits[[1L]], ...)

#' @export
nobs.svyrcs_mifit <- function(object, ...) stats::nobs(object$fits[[1L]], ...)

#' @export
formula.svyrcs_mifit <- function(x, ...) stats::formula(x$fits[[1L]], ...)

#' @export
predict.svyrcs_mifit <- function(object, ...) {
  stop_svyrcs("this is a pooled fit across imputations, not a single model: predicting from it ",
              "would silently use one imputation's data. Use predict() on the svyrcs object ",
              "instead, or reach a component fit through fit$model$fits[[i]].")
}

#' @export
print.svyrcs_mifit <- function(x, ...) {
  cat(sprintf("Pooled fit across %d imputations (Rubin's rules)\n", x$m))
  cat(sprintf("Complete-data degrees of freedom: %s\n", format(x$degf)))
  print(round(stats::coef(x), 5))
  invisible(x)
}

## Barnard-Rubin degrees of freedom for a scalar estimand.
##
## mitools::MIcombine() reports (m - 1)(1 + 1/r)^2, which is unbounded above: against a design with
## 31 degrees of freedom it returned Inf, 756785 and 3457 in testing. Using that for a t quantile
## gives intervals far too narrow. The Barnard-Rubin correction bounds the result by the
## complete-data degrees of freedom, which is what a survey analysis needs.
##
## Vectorised over `u` (within-imputation variance) and `b` (between-imputation variance).
##
## The formula is followed literally, including its behaviour at zero between-imputation variance,
## where it returns (nu_com + 1) nu_com / (nu_com + 3) rather than nu_com itself -- about 6% below
## the complete-data value at nu_com = 31. That is a known conservatism of the approximation, not a
## real loss of information. Earlier versions short-circuited to nu_com there, which made a fit with
## nothing to impute reproduce the complete-data one exactly but silently departed from the
## published rule; matching mice:::barnard.rubin is worth more than that convenience.
barnard_rubin_df <- function(u, b, m, nu_com) {
  n <- max(length(u), length(b))
  u <- rep(u, length.out = n)
  b <- rep(b, length.out = n)
  lambda <- ifelse(u + (1 + 1 / m) * b > 0, (1 + 1 / m) * b / (u + (1 + 1 / m) * b), 0)

  ## With nu_com = Inf there is no complete-data df to bound, so this reduces to Rubin's original
  ## (m - 1)/lambda^2. Guard it explicitly: (Inf + 1)/(Inf + 3) is NaN.
  if (!is.finite(nu_com)) {
    return(ifelse(lambda > 0, (m - 1) / lambda^2, Inf))
  }
  tmp <- (1 - lambda) * (1 + nu_com) * nu_com
  (m - 1) * tmp / ((nu_com + 3) * (m - 1) + lambda^2 * tmp)
}

## D1: the multi-parameter pooling rule of Li, Raghunathan and Rubin (1991), with the small-sample
## denominator degrees of freedom of Reiter (2007).
##
## Two things distinguish it from the obvious construction, and both matter. The statistic uses
## (1 + r) Ubar rather than the pooled total covariance Ubar + (1 + 1/m) B; and the denominator df
## comes from Reiter's formula, which accounts for a finite complete-data df -- exactly the survey
## situation, where the design supplies only a few dozen degrees of freedom.
##
## Earlier versions of this package collapsed r to a scalar and reused the scalar Barnard-Rubin
## formula. That was a plausible-looking construction of my own rather than a named rule, which is
## worse than either: a reader cannot tell what was computed. mitml:::.D1 is the reference
## implementation and the tests check against it.
##
## `r` is the average relative increase in variance due to nonresponse across the block.
d1_test <- function(Qbar, Ubar, B, m, dfcom) {
  k <- length(Qbar)
  Uinv <- tryCatch(solve(Ubar), error = function(e) NULL)
  if (is.null(Uinv)) return(NULL)

  r <- (1 + 1 / m) * sum(diag(B %*% Uinv)) / k
  r <- max(r, 0)
  Ttilde <- (1 + r) * Ubar
  Tinv <- tryCatch(solve(Ttilde), error = function(e) NULL)
  if (is.null(Tinv)) return(NULL)

  stat <- as.numeric(t(Qbar) %*% Tinv %*% Qbar) / k
  t <- k * (m - 1)

  ## Reiter's expansion contains 1/(t - 4) and so is undefined at t = k(m - 1) <= 4 -- reachable
  ## with, say, m = 3 and a two-parameter block. mitml returns NaN there. Fall back to the
  ## complete-data degrees of freedom, which is defined, conservative, and what the reference
  ## distribution tends to when there is little missing information.
  if (is.finite(dfcom) && t <= 4) {
    warning("the small-sample correction for a multi-parameter imputation test needs ",
            "k(m - 1) > 4, but k = ", k, " and m = ", m, " give ", t, ". Using the complete-data ",
            "degrees of freedom instead; increase the number of imputations for a better ",
            "approximation.", call. = FALSE)
    return(list(F = stat, k = k, v = dfcom, r = r))
  }

  v <- if (is.finite(dfcom)) {
    ## Reiter (2007), as implemented in mitml:::.D1.
    a <- r * t / (t - 2)
    vstar <- ((dfcom + 1) / (dfcom + 3)) * dfcom
    c0 <- 1 / (t - 4)
    c1 <- vstar - 2 * (1 + a)
    c2 <- vstar - 4 * (1 + a)
    z <- 1 / c2 +
      c0 * (a^2 * c1 / ((1 + a)^2 * c2)) +
      c0 * (8 * a^2 * c1 / ((1 + a) * c2^2) + 4 * a^2 / ((1 + a) * c2)) +
      c0 * (4 * a^2 / (c2 * c1) + 16 * a^2 * c1 / c2^3) +
      c0 * (8 * a^2 / c2^2)
    4 + 1 / z
  } else if (r <= 0) {
    ## No between-imputation variance and no complete-data df to bound: the chi-square limit.
    Inf
  } else if (t > 4) {
    4 + (t - 4) * (1 + (1 - 2 / t) / r)^2
  } else {
    t * (1 + 1 / k) * (1 + 1 / r)^2 / 2
  }

  list(F = stat, k = k, v = v, r = r)
}

## Fraction of missing information for a scalar estimand.
fraction_missing <- function(u, b, m) {
  total <- u + (1 + 1 / m) * b
  ifelse(total > 0, (1 + 1 / m) * b / total, 0)
}

## The `mi` bundle threaded into the contrast engine and the tests.
## `degf` overrides the complete-data degrees of freedom stored on the pooled fit, so that a `df`
## argument reaches the imputation arithmetic. Without it `df` silently affected only the
## single-design path.
mi_context <- function(fit, degf = NULL) {
  if (!inherits(fit, "svyrcs_mifit")) return(NULL)
  list(Ubar = fit$Ubar, B = fit$B, m = fit$m, degf = degf %||% fit$degf)
}
