model_description <- function(object) {
  if (inherits(object$model, "coxph")) return("survey-weighted Cox proportional hazards")
  fam <- tryCatch(stats::family(object$model), error = function(e) NULL)
  if (is.null(fam)) return("survey-weighted generalised linear model")
  paste0("survey-weighted GLM (", fam$family, ", ", fam$link, " link)")
}

## Contiguous stretches of the curve whose confidence band excludes the null.
##
## Reported because it is the question a reader of the plot actually asks -- "over what range of the
## exposure is the association distinguishable from no effect?" -- and reading it off a figure is
## guesswork.
significant_ranges <- function(curve, null) {
  above <- curve$conf.low > null
  below <- curve$conf.high < null
  out <- list()
  for (dir in c("above", "below")) {
    flag <- if (dir == "above") above else below
    flag[is.na(flag)] <- FALSE
    if (!any(flag)) next
    r <- rle(flag)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    for (i in which(r$values)) {
      out[[length(out) + 1L]] <- data.frame(
        from = curve$x[starts[i]], to = curve$x[ends[i]], direction = dir
      )
    }
  }
  if (!length(out)) {
    return(data.frame(from = numeric(0), to = numeric(0), direction = character(0)))
  }
  res <- do.call(rbind, out)
  res[order(res$from), , drop = FALSE]
}

curve_features <- function(object) {
  cv <- object$curve
  imin <- which.min(cv$estimate)
  imax <- which.max(cv$estimate)
  list(
    min = cv[imin, , drop = FALSE],
    max = cv[imax, , drop = FALSE],
    significant = significant_ranges(cv, object$null)
  )
}

fmt_effect <- function(row, measure, digits = 3) {
  sprintf("%s = %s (%s, %s)", measure,
          fmt_num(row$estimate, digits), fmt_num(row$conf.low, digits),
          fmt_num(row$conf.high, digits))
}

fmt_test <- function(tt, label, width = 20L) {
  if (is.na(tt$F)) {
    return(sprintf("  %-*s not estimable (a linear term only)", width, label))
  }
  sprintf("  %-*s F = %s on %d and %d df,  p %s%s", width, label,
          fmt_num(tt$F, 4), tt$df1, round(tt$df2),
          if (tt$p_F < 2e-16) "" else "= ", fmt_p(tt$p_F))
}

#' @export
print.svyrcs <- function(x, ...) {
  cat("Restricted cubic spline on a complex survey design\n\n")
  cat(sprintf("  Model      %s\n", model_description(x)))
  cat(sprintf("  Outcome    %s\n", deparse1(x$formula[[2L]])))
  cat(sprintf("  Exposure   %s, %d knots at %s%s\n", x$var, x$nk,
              paste(fmt_num(x$knots, 4), collapse = ", "),
              if (isTRUE(x$weighted_knots)) " (weighted quantiles)" else ""))
  cat(sprintf("  Sample     %s observations%s, %d design df\n",
              format(x$n, big.mark = ","),
              if (!is.na(x$nevents)) paste0(", ", format(x$nevents, big.mark = ","), " events") else "",
              round(x$degf)))
  cat(sprintf("  Reference  %s = %s (%s), %s = %s\n", x$var, fmt_num(x$ref$value, 4),
              x$ref$method, x$measure, format(x$null)))
  if (isTRUE(x$n_dropped > 0L)) {
    cat(sprintf("  Dropped    %s rows with missing values\n",
                format(x$n_dropped, big.mark = ",")))
  }
  cat("\n")
  cat(fmt_test(x$tests$overall, "Overall association"), "\n")
  cat(fmt_test(x$tests$nonlinear, "Non-linearity"), "\n")
  invisible(x)
}

#' Summarise a survey restricted cubic spline fit
#'
#' In addition to what [print()] shows, the summary reports the shape of the fitted curve: where the
#' estimated effect is smallest and largest, and over which stretches of the exposure the confidence
#' band excludes the null.
#'
#' @param object An object from [svyrcs()].
#' @param ... Ignored.
#'
#' @return An object of class `summary.svyrcs`, a list with the fit metadata, the two tests and a
#'   `features` component holding the minimum, the maximum and the significant ranges.
#'
#' @examples
#' design <- survey::svydesign(
#'   id = ~psu, strata = ~strata, weights = ~weight,
#'   nest = TRUE, data = nhanes_bmi
#' )
#' fit <- svyrcs(tchol ~ rcs(bmi, 4) + age + sex, design = design)
#' summary(fit)
#'
#' @export
summary.svyrcs <- function(object, ...) {
  structure(
    list(fit = object, features = curve_features(object)),
    class = "summary.svyrcs"
  )
}

#' @export
print.summary.svyrcs <- function(x, ...) {
  fit <- x$fit
  print(fit)
  f <- x$features
  cat("\n  Curve shape\n")
  cat(sprintf("    Lowest    %s = %s,  %s\n", fit$var, fmt_num(f$min$x, 4),
              fmt_effect(f$min, fit$measure)))
  cat(sprintf("    Highest   %s = %s,  %s\n", fit$var, fmt_num(f$max$x, 4),
              fmt_effect(f$max, fit$measure)))

  sig <- f$significant
  if (!nrow(sig)) {
    cat(sprintf("    The %g%% band includes %s = %s across the whole range\n",
                100 * fit$level, fit$measure, format(fit$null)))
  } else {
    cat(sprintf("    %g%% band excludes %s = %s over:\n", 100 * fit$level, fit$measure,
                format(fit$null)))
    for (i in seq_len(nrow(sig))) {
      cat(sprintf("      %s %s to %s  (%s %s %s)\n", fit$var,
                  fmt_num(sig$from[i], 4), fmt_num(sig$to[i], 4), fit$measure,
                  if (sig$direction[i] == "above") ">" else "<", format(fit$null)))
    }
  }
  invisible(x)
}

#' Effect estimates at specific exposure values
#'
#' Gives the estimated effect, with confidence interval, at exposure values you choose, rather than
#' on the plotting grid. This is what you want when a paper reports, say, the hazard ratio at a BMI
#' of 30 against a BMI of 25.
#'
#' @param object An object from [svyrcs()].
#' @param x Numeric vector of exposure values.
#' @param ref Reference value to contrast against. Defaults to the reference the fit was anchored to;
#'   supply a number to compare against a different value without refitting.
#' @param level Confidence level. Defaults to the level used by the fit.
#' @param ... Ignored.
#'
#' @return A data frame with columns `x`, `estimate`, `conf.low`, `conf.high` and `se`.
#'
#' @examples
#' design <- survey::svydesign(
#'   id = ~psu, strata = ~strata, weights = ~weight,
#'   nest = TRUE, data = nhanes_bmi
#' )
#' fit <- svyrcs(tchol ~ rcs(bmi, 4) + age + sex, design = design)
#'
#' # against the fitted reference
#' predict(fit, x = c(20, 25, 30, 35))
#'
#' # against a BMI of 25 instead
#' predict(fit, x = c(20, 30, 35), ref = 25)
#'
#' @export
predict.svyrcs <- function(object, x, ref = NULL, level = NULL, ...) {
  if (missing(x) || !is.numeric(x) || !length(x)) {
    stop_svyrcs("`x` must be a non-empty numeric vector of exposure values")
  }
  svyrcs_curve(
    object$model,
    var = object$var,
    ref = ref %||% object$ref$value,
    at = x,
    level = level %||% object$level,
    design = fit_design(object$model),
    degf = object$degf
  )
}

#' @export
as.data.frame.svyrcs_curve <- function(x, ...) {
  attrs <- setdiff(names(attributes(x)), c("names", "row.names", "class"))
  for (a in attrs) attr(x, a) <- NULL
  class(x) <- "data.frame"
  x
}

#' @export
coef.svyrcs <- function(object, ...) stats::coef(object$model, ...)

#' @export
vcov.svyrcs <- function(object, ...) stats::vcov(object$model, ...)

#' @export
nobs.svyrcs <- function(object, ...) object$n

#' @export
formula.svyrcs <- function(x, ...) x$formula
