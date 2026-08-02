## Harrell's recommended knot placement (Regression Modeling Strategies, 2nd ed., Table 2.3).
harrell_knot_probs <- function(nk) {
  switch(
    as.character(nk),
    "3" = c(0.10, 0.50, 0.90),
    "4" = c(0.05, 0.35, 0.65, 0.95),
    "5" = c(0.05, 0.275, 0.50, 0.725, 0.95),
    "6" = c(0.05, 0.23, 0.41, 0.59, 0.77, 0.95),
    "7" = c(0.025, 0.1833, 0.3417, 0.50, 0.6583, 0.8167, 0.975),
    stop_svyrcs("`knots` must be a count between 3 and 7, or an explicit vector of at least 3 ",
                "knot locations; got ", nk)
  )
}

#' Knot locations for a restricted cubic spline
#'
#' Resolves the `knots` argument of [rcs()] into an explicit, sorted vector of knot locations.
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
  } else {
    if (!is.numeric(knots) || anyNA(knots)) {
      stop_svyrcs("`knots` must be a numeric vector without missing values")
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
  kn
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
#' is Harrell's truncated power parameterisation and is numerically identical to `rms::rcs()`, but it
#' stores its knots on the returned object and registers a [stats::makepredictcall()] method, so a
#' model fitted with `rcs()` reuses exactly the same knots when predicting on new data.
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
#' @section Masking `rms`:
#' Attaching `svyrcs` after `rms` masks `rms::rcs()`. The two bases agree to about 1e-14, so results
#' are unchanged; use `rms::rcs()` explicitly if you need the `rms` version.
#'
#' @seealso [svyrcs()], [rcs_knots()]
#'
#' @examples
#' b <- rcs(nhanes_bmi$bmi, 4)
#' dim(b)
#' attr(b, "knots")
#'
#' @export
rcs <- function(x, knots = 4) {
  var <- deparse1(substitute(x))
  kn <- rcs_knots(x, knots, var = var)
  B <- rcs_design_matrix(as.numeric(x), kn)
  structure(B, knots = kn, nk = length(kn), var = var,
            class = c("svyrcs_basis", "basis", "matrix", "array"))
}

#' Reuse fitted knots when predicting on new data
#'
#' [stats::makepredictcall()] method for the basis returned by [rcs()]. It rewrites the `rcs()` call
#' stored in a model's terms so that it carries the knot locations chosen at fit time. Without it,
#' `predict(fit, newdata = ...)` would re-derive knots from the quantiles of `newdata`, silently
#' evaluating a different spline.
#'
#' @param var The fitted basis, as produced by [rcs()].
#' @param call The call to `rcs()` recorded in the model's terms.
#'
#' @return The call with the fitted knots substituted for whatever `knots` was originally given.
#'
#' @examples
#' fit <- lm(bmi ~ rcs(age, 4), data = nhanes_bmi)
#' # the terms now carry explicit knots rather than "4"
#' attr(terms(fit), "predvars")
#'
#' @export
makepredictcall.svyrcs_basis <- function(var, call) {
  if (as.character(call)[1L] != "rcs") return(call)
  ## Drop every argument except x, so a positional `knots` cannot collide with the one we add back.
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
