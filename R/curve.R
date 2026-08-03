## Locate the rcs() term inside a fitted model.
##
## The knots are read back out of the terms' `predvars`, which makepredictcall.svyrcs_basis() has
## rewritten to carry explicit locations. That makes this work for any model fitted with rcs() in the
## formula, whether or not it went through svyrcs().
find_rcs_term <- function(fit, var = NULL) {
  tt <- stats::terms(fit)
  vars <- attr(tt, "variables")
  pv <- attr(tt, "predvars")
  if (is.null(pv)) pv <- vars

  hits <- list()
  for (i in seq_along(vars)[-1L]) {
    v <- vars[[i]]
    if (!is.call(v) || !identical(as.character(v[[1L]]), "rcs")) next
    pvi <- pv[[i]]
    knots <- tryCatch(eval(pvi$knots), error = function(e) NULL)
    if (is.null(knots) || length(knots) < 3L) {
      stop_svyrcs("the fitted model's rcs() term does not carry explicit knots; refit with ",
                  "svyrcs() or with svyrcs::rcs() so that the knots are recorded")
    }
    hits[[length(hits) + 1L]] <- list(
      label = deparse1(v),
      var = deparse1(v[[2L]]),
      knots = as.numeric(knots),
      nk = length(knots)
    )
  }

  if (!length(hits)) {
    stop_svyrcs("no rcs() term found in the model formula; the exposure must enter the model as, ",
                "for example, rcs(bmi, 4)")
  }
  if (!is.null(var)) {
    keep <- vapply(hits, function(h) identical(h$var, var), logical(1L))
    if (!any(keep)) {
      stop_svyrcs("no rcs() term for '", var, "' in the model; found: ",
                  paste(vapply(hits, `[[`, character(1L), "var"), collapse = ", "))
    }
    hits <- hits[keep]
  }
  if (length(hits) > 1L) {
    stop_svyrcs("the model has ", length(hits), " rcs() terms (",
                paste(vapply(hits, `[[`, character(1L), "var"), collapse = ", "),
                "); pass `var` to choose one")
  }
  hits[[1L]]
}

## Column names a k-knot basis contributes: the term label, then the label with primes.
spline_colnames <- function(label, nk) {
  paste0(label, c("", vapply(seq_len(nk - 2L), function(i) strrep("'", i), character(1L))))
}

## Exact-match the spline coefficients.
##
## Matching by name rather than by grepping the variable name matters: a covariate called
## `bmi_group` would be swept into a grep for "bmi", silently corrupting every estimate.
spline_coef_index <- function(fit, term) {
  want <- spline_colnames(term$label, term$nk)
  idx <- match(want, names(stats::coef(fit)))
  if (anyNA(idx)) {
    stop_svyrcs("could not find the spline coefficients for '", term$var, "' in the fitted model; ",
                "expected ", paste(sQuote(want[is.na(idx)]), collapse = ", "), ". If the spline is ",
                "interacted with a modifier, write rcs(", term$var, ", k) * m so that the main ",
                "effects are estimated too.")
  }
  idx
}

## The survey design a model was fitted to, if it kept one.
fit_design <- function(fit) {
  d <- fit$survey.design
  if (is.null(d)) d <- fit$design
  if (inherits(d, "survey.design") || inherits(d, "svyrep.design")) d else NULL
}

## Effect measure and whether the linear-predictor contrast should be exponentiated.
effect_measure <- function(fit) {
  if (inherits(fit, "coxph")) {
    return(list(measure = "HR", exponentiate = TRUE, null = 1))
  }
  fam <- tryCatch(stats::family(fit), error = function(e) NULL)
  link <- if (is.null(fam)) "identity" else fam$link
  fname <- if (is.null(fam)) "" else fam$family
  switch(
    link,
    logit = list(measure = "OR", exponentiate = TRUE, null = 1),
    ## A log link on a count or binary family gives a risk or rate ratio; on a gaussian one it is a
    ## ratio of means, which is not the same thing and should not be labelled RR.
    log = if (grepl("poisson|binomial", fname)) {
      list(measure = "RR", exponentiate = TRUE, null = 1)
    } else {
      list(measure = "Ratio", exponentiate = TRUE, null = 1)
    },
    identity = list(measure = "Difference", exponentiate = FALSE, null = 0),
    {
      warning("no standard effect measure for the '", link, "' link; reporting contrasts on the ",
              "link scale", call. = FALSE)
      list(measure = paste0(link, "-scale difference"), exponentiate = FALSE, null = 0)
    }
  )
}

## The numerical core.
##
## Because the linear predictor is additive and every covariate takes the same value at x and at the
## reference x0, all non-spline columns cancel in the difference. Only the spline block is needed:
##
##   dB     = B(x) - B(x0)
##   effect = dB %*% beta,  var = rowSums((dB %*% V) * dB)
##
## This is exact, not an approximation: it reproduces predict(type = "lp") differences to ~1e-15.
contrast_estimates <- function(x, x0, knots, beta, V, degf, level, exponentiate, L = NULL,
                               mi = NULL) {
  B <- rcs_design_matrix(x, knots)
  B0 <- rcs_design_matrix(x0, knots)
  dB <- sweep(B, 2L, as.numeric(B0), "-")

  ## `L` maps the full coefficient vector to the spline coefficients in force for this group: the
  ## main effects alone at the reference level, main + interaction elsewhere. Doing it as one matrix
  ## rather than summing coefficient vectors is what makes the variance pick up
  ## Cov(main, interaction); adding separate variances would understate every non-reference SE.
  M <- if (is.null(L)) dB else dB %*% L

  est <- as.numeric(M %*% beta)
  variance <- rowSums((M %*% V) * M)
  ## tiny negative values can appear at the reference point through cancellation
  variance[variance < 0 & variance > -1e-10] <- 0
  se <- sqrt(variance)

  ## Under multiple imputation the reference distribution is not the same at every point: the
  ## between-imputation variance of the contrast varies along the curve, so the degrees of freedom
  ## have to be computed per point rather than once for the whole fit.
  extra <- NULL
  if (is.null(mi)) {
    df_point <- degf
  } else {
    u <- rowSums((M %*% mi$Ubar) * M)
    b <- rowSums((M %*% mi$B) * M)
    u[u < 0] <- 0
    b[b < 0] <- 0
    df_point <- barnard_rubin_df(u, b, mi$m, mi$degf)
    extra <- list(df = df_point, fmi = fraction_missing(u, b, mi$m))
  }

  crit <- stats::qt(1 - (1 - level) / 2, df = df_point)
  lo <- est - crit * se
  hi <- est + crit * se

  out <- if (exponentiate) {
    data.frame(x = x, estimate = exp(est), conf.low = exp(lo), conf.high = exp(hi),
               se = se, row.names = NULL)
  } else {
    data.frame(x = x, estimate = est, conf.low = lo, conf.high = hi,
               se = se, row.names = NULL)
  }
  if (!is.null(extra)) {
    out$df <- extra$df
    out$fmi <- extra$fmi
  }
  out
}

#' Exposure-response curve from a fitted survey model
#'
#' Estimates the exposure-response curve implied by an `rcs()` term in an already-fitted survey
#' model, as a contrast against a reference exposure value. [svyrcs()] calls this internally; use it
#' directly when you have fitted the model yourself, for instance because you need a model
#' specification `svyrcs()` does not build for you.
#'
#' The curve is computed from the spline coefficients alone. Every covariate takes the same value at
#' `x` and at the reference, so covariate contributions cancel exactly and no covariate reference
#' grid is needed. Confidence intervals use a *t* quantile on the survey degrees of freedom.
#'
#' @param fit A fitted model containing an [rcs()] term, typically from [survey::svycoxph()] or
#'   [survey::svyglm()]. Any model with `coef()`, `vcov()` and `terms()` methods works.
#' @param var Name of the exposure. Only needed when the model contains more than one `rcs()` term.
#' @param ref Reference value: a number, or one of `"median"`, `"mean"`, `"quantile"`, `"min"`
#'   (minimum-risk point) or `"max"` (maximum-risk point). See [svyrcs()].
#' @param ref_prob Probability used when `ref = "quantile"`.
#' @param at Optional numeric vector of exposure values to evaluate. When given, `n` and `range` are
#'   ignored and the result has one row per element of `at`.
#' @param n Number of points on the curve.
#' @param range Length-2 numeric giving the exposure range to plot over. Defaults to the 1st and 99th
#'   survey-weighted percentiles of the exposure.
#' @param level Confidence level.
#' @param design The survey design. Defaults to the design stored on `fit`, which is what
#'   [survey::svyglm()] and [survey::svycoxph()] keep. Needed for weighted reference values and
#'   ranges.
#' @param degf Design degrees of freedom for the *t* quantile. Defaults to `survey::degf(design)`.
#' @param group When the spline is interacted with an effect modifier, the level or levels to
#'   estimate. `NULL` (default) returns every level, stacked, with a `group` column.
#'
#' @return A data frame of class `svyrcs_curve` with columns `x`, `estimate`, `conf.low`,
#'   `conf.high` and `se` (on the linear-predictor scale), plus `group` when the model has an
#'   effect modifier. Carries attributes `ref`, `ref_method`, `measure`, `null`, `knots`, `var`,
#'   `degf`, `level` and `modifier`.
#'
#' @seealso [svyrcs()]
#'
#' @examples
#' design <- survey::svydesign(
#'   id = ~psu, strata = ~strata, weights = ~weight,
#'   nest = TRUE, data = nhanes_bmi
#' )
#' fit <- survey::svyglm(tchol ~ rcs(bmi, 4) + age + sex, design = design)
#' curve <- svyrcs_curve(fit, "bmi", n = 20)
#' head(curve)
#'
#' @export
svyrcs_curve <- function(fit, var = NULL, ref = "median", ref_prob = 0.5, at = NULL, n = 200,
                         range = NULL, level = 0.95, design = NULL, degf = NULL, group = NULL) {
  term <- find_rcs_term(fit, var)
  spline_coef_index(fit, term)  # validates that the main effects are present
  modifier <- find_modifier(fit, term)
  ## The full coefficient vector, with the selection matrix doing the extraction: one code path for
  ## grouped and ungrouped fits.
  beta <- stats::coef(fit)
  V <- stats::vcov(fit)
  meas <- effect_measure(fit)

  design <- design %||% fit_design(fit)
  ## Weighted quantities are averaged over the imputations, so the design is always handled as a
  ## list. For an ordinary fit that list has one element and nothing changes.
  designs <- as_design_list(design)
  degf <- degf %||% (if (inherits(fit, "svyrcs_mifit")) fit$degf
                     else if (!is.null(designs)) survey::degf(designs[[1L]]) else NULL)
  if (is.null(degf) || !is.finite(degf) || degf <= 0) {
    stop_svyrcs("could not determine the survey degrees of freedom; pass `degf` explicitly")
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop_svyrcs("`level` must be a single number strictly between 0 and 1")
  }

  xvals <- exposure_values(fit, term$var, designs)
  rng <- range %||% exposure_range(xvals, designs, term$var)
  grid <- if (!is.null(at)) {
    if (!is.numeric(at) || !length(at)) stop_svyrcs("`at` must be a non-empty numeric vector")
    as.numeric(at)
  } else {
    if (!is_count(n) || n < 2) stop_svyrcs("`n` must be a whole number of at least 2; got ", n)
    seq(rng[1L], rng[2L], length.out = as.integer(n))
  }

  ## A restricted cubic spline is linear beyond the outer knots, so extrapolation is defined -- which
  ## is exactly what makes it easy to do by accident. `ref` is refused outright, because anchoring
  ## outside the data makes every estimate an extrapolation; individual points only warn, because
  ## deliberate extrapolation is a legitimate, if bold, choice.
  warn_outside_range(grid, xvals, term$var, if (is.null(at)) "`range`" else "`at`")

  wanted <- resolve_groups(modifier, group)

  ## Flag a design that cannot identify the spline, here as well as in the tests, so the warning
  ## does not depend on which entry point the user came through.
  sp_idx <- match(spline_colnames(term$label, term$nk), names(beta))
  warn_if_rank_deficient(V[sp_idx, sp_idx, drop = FALSE], term$var)

  mi <- mi_context(fit)
  curve_at <- function(x0, xs, g = NULL) {
    L <- group_selection(names(beta), term, modifier, g)
    contrast_estimates(xs, x0, term$knots, beta, V, degf, level, meas$exponentiate, L, mi)
  }

  ## One reference for all groups: curves anchored at different exposures are not comparable, and
  ## comparability is the entire point of estimating them together. A curve-derived reference
  ## ("min" / "max") is therefore located on the reference level's curve.
  ref_group <- if (is.null(modifier)) NULL else modifier$ref_level
  refinfo <- resolve_ref(ref, xvals, designs, term$var, ref_prob = ref_prob,
                         curve_fun = function(x0, xs) curve_at(x0, xs, ref_group),
                         range = rng, meas = meas)
  if (!is.null(modifier) && refinfo$method %in% c("minimum-risk point", "maximum-risk point")) {
    refinfo$method <- paste0(refinfo$method, " (", modifier$ref_level, ")")
  }

  out <- if (is.null(modifier)) {
    curve_at(refinfo$value, grid)
  } else {
    blocks <- lapply(wanted, function(g) {
      cbind(curve_at(refinfo$value, grid, g), group = g, stringsAsFactors = FALSE)
    })
    res <- do.call(rbind, blocks)
    res$group <- factor(res$group, levels = modifier$levels)
    res
  }

  if (any(!is.finite(out$estimate))) {
    warning("the curve contains non-finite estimates, which means the underlying model did not ",
            "converge to finite coefficients -- perfect separation is the usual cause. Check the ",
            "fitted model before using these results.", call. = FALSE)
  }

  structure(
    out,
    ref = refinfo$value, ref_method = refinfo$method,
    measure = meas$measure, null = meas$null, exponentiate = meas$exponentiate,
    knots = term$knots, var = term$var, label = term$label, degf = degf, level = level,
    modifier = modifier,
    class = c("svyrcs_curve", "data.frame")
  )
}

## Which modifier levels to estimate, validating anything the caller asked for.
resolve_groups <- function(modifier, group) {
  if (is.null(modifier)) {
    if (!is.null(group)) {
      stop_svyrcs("`group` was given but the model has no effect modifier interacting with the ",
                  "spline term")
    }
    return(NULL)
  }
  if (is.null(group)) return(modifier$levels)
  group <- as.character(group)
  bad <- setdiff(group, modifier$levels)
  if (length(bad)) {
    stop_svyrcs("unknown level", if (length(bad) > 1L) "s" else "", " of '", modifier$var, "': ",
                paste(sQuote(bad), collapse = ", "), ". Available: ",
                paste(sQuote(modifier$levels), collapse = ", "))
  }
  group
}

warn_outside_range <- function(x, xvals, var, what) {
  obs <- range(xvals, na.rm = TRUE)
  out <- x < obs[1L] | x > obs[2L]
  if (!any(out, na.rm = TRUE)) return(invisible(FALSE))
  warning(sum(out, na.rm = TRUE), " of ", length(x), " values in ", what, " fall outside the ",
          "observed range of '", var, "' (", format(obs[1L]), " to ", format(obs[2L]), "). ",
          "The spline is linear beyond the outer knots, so these are extrapolations.",
          call. = FALSE)
  invisible(TRUE)
}

## Normalise a design argument to a list of component designs. An ordinary design becomes a
## one-element list; a multiply imputed one becomes its components.
as_design_list <- function(design) {
  if (is.null(design)) return(NULL)
  if (is_mi_design(design)) return(mi_components(design))
  if (is.list(design) && !inherits(design, c("survey.design", "svyrep.design")) &&
      length(design) && inherits(design[[1L]], c("survey.design", "svyrep.design"))) {
    return(design)
  }
  list(design)
}

## Exposure values actually used in the fit, for range and reference calculations. Under imputation
## the exposure itself may be imputed, so every completed version contributes.
exposure_values <- function(fit, var, designs) {
  if (!is.null(designs)) {
    v <- tryCatch(unlist(lapply(designs, function(d) d$variables[[var]]), use.names = FALSE),
                  error = function(e) NULL)
    if (!is.null(v) && length(v)) return(as.numeric(v))
  }
  mf <- tryCatch(stats::model.frame(fit), error = function(e) NULL)
  if (!is.null(mf)) {
    for (nm in names(mf)) {
      col <- mf[[nm]]
      if (inherits(col, "svyrcs_basis") && identical(attr(col, "var"), var)) {
        return(as.numeric(col[, 1L]))
      }
      if (identical(nm, var) && is.numeric(col)) return(as.numeric(col))
    }
  }
  stop_svyrcs("could not recover the values of '", var, "' from the model or design; ",
              "pass `range` explicitly")
}

exposure_range <- function(xvals, designs, var) {
  if (!is.null(designs)) {
    q <- tryCatch(pooled_weighted_quantile(var, designs, c(0.01, 0.99)), error = function(e) NULL)
    if (!is.null(q)) return(unname(q))
  }
  unname(stats::quantile(xvals, c(0.01, 0.99), na.rm = TRUE))
}
