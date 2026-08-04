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
## Vectorised over `u` (within-imputation variance) and `b` (between-imputation variance). Returns
## nu_com exactly where b <= 0, so a fit with nothing to impute is exact rather than merely close.
barnard_rubin_df <- function(u, b, m, nu_com) {
  out <- rep(nu_com, length.out = max(length(u), length(b)))
  live <- b > 0 & u > 0
  ## With nu_com = Inf the complete-data correction has nothing to bound, so Barnard-Rubin reduces to
  ## Rubin's original (m - 1)(1 + 1/r)^2. Guard it explicitly: (Inf + 1)/(Inf + 3) is NaN.
  if (!is.finite(nu_com)) {
    if (!any(live)) return(out)
    r <- (1 + 1 / m) * b[live] / u[live]
    out[live] <- (m - 1) * (1 + 1 / r)^2
    return(out)
  }
  if (!any(live)) return(out)

  bb <- b[live]
  uu <- u[live]
  r <- (1 + 1 / m) * bb / uu
  nu_old <- (m - 1) * (1 + 1 / r)^2
  gamma <- (1 + 1 / m) * bb / (uu + (1 + 1 / m) * bb)
  nu_obs <- (nu_com + 1) / (nu_com + 3) * nu_com * (1 - gamma)
  out[live] <- nu_old * nu_obs / (nu_old + nu_obs)
  out
}

## Fraction of missing information for a scalar estimand.
fraction_missing <- function(u, b, m) {
  total <- u + (1 + 1 / m) * b
  ifelse(total > 0, (1 + 1 / m) * b / total, 0)
}

## The `mi` bundle threaded into the contrast engine and the tests.
mi_context <- function(fit) {
  if (!inherits(fit, "svyrcs_mifit")) return(NULL)
  list(Ubar = fit$Ubar, B = fit$B, m = fit$m, degf = fit$degf)
}
