#' Reference values for an exposure-response curve
#'
#' A restricted cubic spline curve is only identified up to a reference value: every estimate is a
#' contrast against some exposure level, which is where the curve passes through the null. `svyrcs`
#' supports the reference choices commonly seen in the survey literature.
#'
#' @section Available methods:
#' \describe{
#'   \item{a number}{Used as given. Must lie inside the observed exposure range.}
#'   \item{`"median"`}{The survey-weighted median of the exposure. The default.}
#'   \item{`"mean"`}{The survey-weighted mean.}
#'   \item{`"quantile"`}{The survey-weighted quantile at `ref_prob`.}
#'   \item{`"min"`}{The exposure value at which the fitted curve is lowest: the minimum-risk point.}
#'   \item{`"max"`}{The exposure value at which the fitted curve is highest.}
#' }
#'
#' Weighted quantiles and means are computed from the survey design with [survey::svyquantile()] and
#' [survey::svymean()], so the reference is a population quantity rather than a property of the
#' unweighted sample. `"min"` and `"max"` are found by evaluating the curve on a dense grid against a
#' provisional (weighted median) reference; because the reference only shifts the curve, the located
#' point does not depend on that provisional choice.
#'
#' @name references
NULL

## Returns list(value, method). `curve_fun(x0, xs)` must return a data frame with columns x and
## estimate, i.e. the curve referenced at x0.
## `designs` is always a list of component designs -- length 1 for an ordinary fit, m for a multiply
## imputed one. Weighted references are averaged over it so that every imputation is anchored at the
## same value.
resolve_ref <- function(ref, xvals, designs, var, ref_prob = 0.5, curve_fun = NULL,
                        range = NULL, meas = NULL, grid_n = 500L) {
  xr <- range(xvals, na.rm = TRUE)

  if (is.numeric(ref)) {
    if (length(ref) != 1L || is.na(ref)) {
      stop_svyrcs("a numeric `ref` must be a single non-missing value")
    }
    if (ref < xr[1L] || ref > xr[2L]) {
      stop_svyrcs("`ref` = ", format(ref), " is outside the observed range of '", var, "' (",
                  format(xr[1L]), " to ", format(xr[2L]), ")")
    }
    return(list(value = as.numeric(ref), method = "user-specified"))
  }

  if (!is.character(ref) || length(ref) != 1L) {
    stop_svyrcs("`ref` must be a single number, or one of \"median\", \"mean\", \"quantile\", ",
                "\"min\" or \"max\"")
  }

  allowed <- c("median", "mean", "quantile", "min", "max")
  if (!ref %in% allowed) {
    stop_svyrcs("`ref` must be a number, or one of ", paste(dQuote(allowed, FALSE), collapse = ", "),
                "; got ", dQuote(ref, FALSE))
  }

  if (ref %in% c("median", "quantile")) {
    p <- if (ref == "median") 0.5 else ref_prob
    if (!is.numeric(p) || length(p) != 1L || is.na(p) || p <= 0 || p >= 1) {
      stop_svyrcs("`ref_prob` must be a single number strictly between 0 and 1; got ", format(p))
    }
    value <- if (!is.null(designs)) {
      unname(pooled_weighted_quantile(var, designs, p))
    } else {
      unname(stats::quantile(xvals, p, na.rm = TRUE))
    }
    method <- if (ref == "median") "weighted median" else paste0("weighted ", p * 100, "th percentile")
    if (is.null(designs)) method <- sub("weighted", "unweighted", method)
    return(list(value = value, method = method))
  }

  if (ref == "mean") {
    value <- if (!is.null(designs)) pooled_weighted_mean(var, designs) else mean(xvals, na.rm = TRUE)
    return(list(value = value,
                method = if (is.null(designs)) "unweighted mean" else "weighted mean"))
  }

  ## "min" / "max": locate the extremum of the curve itself.
  if (is.null(curve_fun)) {
    stop_svyrcs("`ref = \"", ref, "\"` needs the fitted curve; use svyrcs() or svyrcs_curve()")
  }
  rng <- range %||% unname(stats::quantile(xvals, c(0.01, 0.99), na.rm = TRUE))
  provisional_x0 <- if (!is.null(designs)) {
    unname(pooled_weighted_quantile(var, designs, 0.5))
  } else {
    stats::median(xvals, na.rm = TRUE)
  }
  grid <- seq(rng[1L], rng[2L], length.out = grid_n)
  est <- curve_fun(provisional_x0, grid)$estimate
  pick <- if (ref == "min") which.min(est) else which.max(est)
  list(value = grid[pick],
       method = if (ref == "min") "minimum-risk point" else "maximum-risk point")
}
