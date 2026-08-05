## Harrell's recommended knot placement (Regression Modeling Strategies, 2nd ed., Table 2.3).
##
## Computed rather than tabulated. The published table rounds to four decimals, and using those
## rounded values put the 7-knot placement about 0.004 away from `Hmisc::rcspline.eval()`; the table
## is itself `seq(outer, 1 - outer, length.out = nk)`, so deriving it removes the discrepancy.
harrell_knot_probs <- function(nk) {
  if (!nk %in% 3:7) {
    stop_svyrcs("`knots` must be a count between 3 and 7, or an explicit vector of at least 3 ",
                "knot locations; got ", nk)
  }
  outer <- if (nk > 6L) 0.025 else if (nk > 3L) 0.05 else 0.1
  seq(outer, 1 - outer, length.out = nk)
}

#' Knot locations for a restricted cubic spline
#'
#' Resolves the `knots` argument of [rcspline()] into an explicit, sorted vector of knot locations.
#' Exported so that knot placement can be inspected, or computed once and reused across models.
#'
#' @param x Numeric vector of exposure values.
#' @param knots Either a single count between 3 and 7, in which case knots are placed at Harrell's
#'   recommended quantiles of `x`, or a numeric vector of at least three explicit knot locations.
#' @param var Optional variable name, used only to make error messages readable.
#'
#' @return A sorted numeric vector of knot locations.
#'
#' @examples
#' rcs_knots(nhanes_bmi$bmi, 4)
#' rcs_knots(nhanes_bmi$bmi, c(20, 25, 30, 40))
#'
#' @export
rcs_knots <- function(x, knots = 4, var = "x") {
  if (!is.numeric(x)) {
    stop_svyrcs("the exposure '", var, "' must be numeric to fit a restricted cubic spline; got ",
                class(x)[1L])
  }
  xo <- x[!is.na(x)]
  if (!length(xo)) stop_svyrcs("the exposure '", var, "' is entirely missing")

  if (length(knots) == 1L) {
    if (!is.numeric(knots) || !is.finite(knots)) {
      stop_svyrcs("`knots` of length 1 is read as the number of knots and must be a finite whole ",
                  "number between 3 and 7; got ", format(knots))
    }
    if (!is_count(knots)) {
      stop_svyrcs("`knots` of length 1 is read as the number of knots and must be a whole number ",
                  "between 3 and 7; got ", format(knots))
    }
    nk <- as.integer(round(knots))
    ## Only meaningful when knots are derived from the data: explicit knots may legitimately be
    ## evaluated on a handful of new observations (e.g. a prediction grid of three points).
    if (length(unique(xo)) < nk) {
      stop_svyrcs("placing ", nk, " knots needs at least ", nk, " distinct values of '", var,
                  "', but it has ", length(unique(xo)))
    }
    kn <- unname(stats::quantile(xo, harrell_knot_probs(nk), na.rm = TRUE, names = FALSE))
    kn <- small_sample_outer_knots(kn, xo, var = var)
  } else {
    if (!is.numeric(knots) || !all(is.finite(knots))) {
      stop_svyrcs("`knots` must be a numeric vector of finite values; got ",
                  paste(format(knots[!is.finite(knots)]), collapse = ", "))
    }
    kn <- sort(unique(as.numeric(knots)))
    nk <- length(kn)
    if (nk < 3L) {
      stop_svyrcs("at least 3 distinct knots are needed for a restricted cubic spline; got ", nk)
    }
  }

  kn <- sort(kn)
  if (anyDuplicated(kn)) {
    stop_svyrcs("could not place ", nk, " distinct knots for '", var, "': the requested quantiles ",
                "coincide, which happens when the exposure has few distinct values. Use fewer ",
                "knots or supply explicit knot locations.")
  }
  warn_if_knots_crowded(kn, var)
  kn
}

## Knots that are distinct but nearly coincident.
##
## `anyDuplicated()` catches exact ties, which is the common case: an exposure with most of its mass
## at one value cannot place distinct quantiles at all. What it cannot catch is a cluster carrying
## floating-point jitter -- a biomarker below the limit of detection substituted by a constant, say,
## where the substitution arithmetic leaves the values differing around 1e-15. Two quantiles then land
## inside that cluster a few 1e-16 apart, distinct enough to pass.
##
## `rcs_design_matrix()` divides by (t_k - t_1)^2 and forms (t_k - t_j)/(t_k - t_{k-1}), so a gap that
## small produces enormous factors and cancellation. Measured against splines::ns() on identical
## knots, the fitted curve is accurate to 4e-12 at a relative gap of 1e-6 but wrong by 2% at 1e-16.
##
## The rank-deficiency warning already fires from about 1e-6, so nothing here is silent. It says the
## design cannot identify the model and suggests more primary sampling units, which will not help a
## spacing problem. This names the actual cause.
##
## The threshold sits where the degradation starts, not where it becomes material, because by the
## time it is material the numbers are already unusable. It does not fire on ordinary exposures --
## uniform, normal, gamma, exponential, lognormal and the shipped NHANES BMI are all well clear.
warn_if_knots_crowded <- function(kn, var, tol = 1e-6) {
  span <- diff(range(kn))
  if (!is.finite(span) || span <= 0) return(invisible(FALSE))
  gaps <- diff(kn)
  i <- which.min(gaps)
  if (gaps[i] / span >= tol) return(invisible(FALSE))
  warning("two knots for '", var, "' are almost coincident: ", format(kn[i], digits = 15), " and ",
          format(kn[i + 1L], digits = 15), " differ by ", format(gaps[i], digits = 3),
          ", which is ", format(gaps[i] / span, digits = 3), " of the knot range. The spline basis ",
          "divides by these gaps, so the fitted curve loses accuracy. This usually means the ",
          "exposure has a dense cluster -- values below a limit of detection replaced by a constant ",
          "are the common case. Use fewer knots, or supply explicit knot locations.", call. = FALSE)
  invisible(TRUE)
}

## Harrell's small-sample adjustment: below 100 observations the outer knots are moved to the fifth
## smallest and fifth largest values, because an extreme quantile of a small sample is an unstable
## place to anchor a spline tail.
##
## This reproduces `Hmisc::rcspline.eval()`'s `length(xx) < 100` branch, except at n = 9 with 4 or 5
## knots, where the replacements collide with the interior knots: Hmisc silently returns three, we
## error. Interior knots are untouched,
## which is why the two agreed on those all along; only the outer pair moved. Adopted so that the
## default placement matches `rms`, not only the basis given identical knots.
##
## At very small n the replacements can fall inside the interior knots -- for `1:10` the result sorts
## to 4.15, 5, 6, 6.85. `Hmisc` behaves the same way, and diverging here would reintroduce exactly
## the sort of private rule this change exists to remove.
small_sample_outer_knots <- function(kn, xo, var = "x") {
  if (length(xo) >= 100L) return(kn)
  nk <- length(kn)
  ## The rule indexes the 5th and (n-4)th order statistics, so it needs at least 5 observations.
  ## Before 0.6.4 `xs[5L]` on a shorter vector gave NA or a zero-length replacement, and the caller
  ## got malformed knots or a base error about replacement length.
  ##
  ## Not 9: at n = 8 the two indices cross (5 and 4), but `sort()` puts them back in order and the
  ## result is still four distinct knots matching Hmisc. A guard set by when the indices look tidy
  ## rather than by when they are defined would reject placements that work.
  xs <- sort(xo)
  if (length(xs) < 5L) {
    stop_svyrcs("placing ", nk, " knots by Harrell's small-sample rule needs at least 5 ",
                "observations of '", var, "', but there are ", length(xs),
                ". Supply explicit knot locations, or use more data.")
  }
  kn[1L] <- xs[5L]
  kn[nk] <- xs[length(xs) - 4L]
  out <- sort(unique(kn))
  ## The replacements can collide with an interior knot, and `unique()` then silently shortens the
  ## vector: `rcs_knots(1:9, 5)` used to return three knots, so the caller fitted a different spline
  ## than it asked for with nothing to show it.
  if (length(out) != nk) {
    stop_svyrcs("could not place ", nk, " distinct knots for '", var, "': below 100 observations ",
                "the outer knots move to the fifth smallest and fifth largest values, and here ",
                "they coincide with the interior knots, leaving ", length(out),
                ". Use fewer knots, or supply explicit knot locations.")
  }
  out
}

## The numerical core: Harrell's truncated power basis, normalised by (t_k - t_1)^2.
## This is rms::rcs()'s default (norm = 2) parameterisation; the two agree to ~1e-14.
rcs_design_matrix <- function(x, knots) {
  k <- length(knots)
  t1 <- knots[1L]
  tkm1 <- knots[k - 1L]
  tk <- knots[k]
  denom <- (tk - t1)^2
  pos3 <- function(u) pmax(u, 0)^3

  B <- matrix(NA_real_, nrow = length(x), ncol = k - 1L)
  B[, 1L] <- x
  for (j in seq_len(k - 2L)) {
    tj <- knots[j]
    B[, j + 1L] <- (pos3(x - tj) -
                      pos3(x - tkm1) * (tk - tj) / (tk - tkm1) +
                      pos3(x - tk) * (tkm1 - tj) / (tk - tkm1)) / denom
  }
  colnames(B) <- c("", vapply(seq_len(k - 2L), function(i) strrep("'", i), character(1L)))
  B
}

#' Restricted cubic spline basis
#'
#' Builds a restricted cubic spline (natural spline) basis for use inside a model formula. The basis
#' is Harrell's truncated power parameterisation and matches `rms::rcs()` -- both the basis and, since
#' 0.6.2, the knot placement this function performs, including the small-sample rule that moves the
#' outer knots to
#' the fifth smallest and fifth largest values below 100 observations. Unlike `rms`, it stores its
#' knots on the returned object and registers a [stats::makepredictcall()] method, so a model fitted
#' with `rcspline()` reuses exactly the same knots when predicting on new data.
#'
#' Because a restricted cubic spline with `k` knots contributes `k - 1` columns, the first of which
#' is `x` itself, a test that the remaining columns are jointly zero is a test of non-linearity. This
#' is what [svyrcs()] reports as the non-linearity *p*-value.
#'
#' @param x Numeric exposure vector.
#' @param knots Either the number of knots (3 to 7, default 4), placed at Harrell's recommended
#'   quantiles of `x`, or an explicit numeric vector of knot locations. [svyrcs()] normally replaces
#'   a count with survey-weighted quantiles before fitting.
#'
#' @return A numeric matrix with `length(x)` rows and `k - 1` columns, of class `svyrcs_basis`, with
#'   attributes `knots` and `nk`.
#'
#' @section Relationship to `rms`:
#' This function used to be called `rcs()`, which masked `rms::rcs()` whenever both packages were
#' attached. It was renamed in 0.7.0 because masking could not be made safe: an `ols()`, `lrm()` or
#' `cph()` model written with a bare `rcs()` would pick up this basis, fit correctly, and then lose
#' the `Nonlinear` row from `anova()` with no error and no warning -- and that row is usually the
#' reason for fitting a spline at all.
#'
#' With the rename there is no collision. `rcs()` in an `rms` formula is unambiguously `rms::rcs()`.
#'
#' The basis itself is Harrell's, and agrees with `rms::rcs()` to about 1e-14 **given the same
#' knots**. Note the converse still holds: passing `rcspline()` into an `rms` model fits, but the
#' result is not an `rms` design term, so `anova()` reports no `Nonlinear` row and `Predict()` fails.
#' Use `rms::rcs()` for `rms` models and `rcspline()` for this package.
#'
#' @seealso [svyrcs()], [rcs_knots()], and `rms::rcs()` for `rms` models -- the two are not
#'   interchangeable in an `rms` formula, see the section above.
#'
#' @examples
#' b <- svyrcs::rcspline(nhanes_bmi$bmi, 4)
#' dim(b)
#' attr(b, "knots")
#'
#' @export
rcspline <- function(x, knots = 4) {
  var <- deparse1(substitute(x))
  kn <- rcs_knots(x, knots, var = var)
  B <- rcs_design_matrix(as.numeric(x), kn)
  structure(B, knots = kn, nk = length(kn), var = var,
            class = c("svyrcs_basis", "basis", "matrix", "array"))
}

#' Reuse fitted knots when predicting on new data
#'
#' [stats::makepredictcall()] method for the basis returned by [rcspline()]. It rewrites the `rcspline()` call
#' stored in a model's terms so that it carries the knot locations chosen at fit time. Without it,
#' `predict(fit, newdata = ...)` would re-derive knots from the quantiles of `newdata`, silently
#' evaluating a different spline.
#'
#' @param var The fitted basis, as produced by [rcspline()].
#' @param call The call to `rcspline()` recorded in the model's terms.
#'
#' @return The call with the fitted knots substituted for whatever `knots` was originally given.
#'
#' @examples
#' fit <- lm(bmi ~ svyrcs::rcspline(age, 4), data = nhanes_bmi)
#' # the terms now carry explicit knots rather than "4"
#' attr(terms(fit), "predvars")
#'
#' @export
makepredictcall.svyrcs_basis <- function(var, call) {
  ## Match the namespaced forms too. Matching the bare name only meant `lm(y ~ svyrcs::rcspline(x, 4))`
  ## recorded the count-based call and re-derived knots from newdata: a model trained on an exposure
  ## of 1-100 predicted a grid near 1000 using knots near 1000. svyrcs() itself was unaffected,
  ## because it inlines the knots before fitting -- which is exactly why the tests missed this.
  if (!is_call_to(call, "rcspline")) return(call)
  ## Keep element 1 (the possibly-qualified function) and element 2 (x), dropping everything else so
  ## a positional `knots` cannot collide with the one added back.
  out <- call[1L:2L]
  out$knots <- attr(var, "knots")
  out
}

#' @export
`[.svyrcs_basis` <- function(x, i, j, ..., drop = FALSE) {
  keep <- attributes(x)[c("knots", "nk", "var")]
  out <- NextMethod("[")
  ## Row subsetting (na.omit, subset) must not silently discard the knots; column subsetting or
  ## dropping to a vector means the object is no longer a basis, so let it degrade.
  if (is.matrix(out) && ncol(out) == ncol(x)) {
    attributes(out) <- c(attributes(out), keep)
    class(out) <- class(x)
  }
  out
}
