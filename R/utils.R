`%||%` <- function(x, y) if (is.null(x)) y else x

## Errors raised by this package carry the class "svyrcs_error" so callers (and tests) can
## distinguish a misuse of svyrcs from a failure inside survey / survival.
stop_svyrcs <- function(...) {
  msg <- paste0(...)
  stop(structure(
    class = c("svyrcs_error", "error", "condition"),
    list(message = msg, call = sys.call(-1L))
  ))
}

## Rows a model fitted from `formula` on `data` would actually use.
##
## Evaluating the model frame is the point: it applies the formula's transformations, so a value that
## only becomes non-finite after evaluation is caught. `complete.cases()` on the raw columns is not
## equivalent and was the source of a silent knot-placement bug.
##
## `na.pass` keeps every row so the result lines up with `data`; non-finite numerics are then treated
## as unusable, which `is.na()` alone would miss for NaN produced by, say, log() of a negative value.
analytic_rows <- function(formula, data) {
  mf <- tryCatch(
    stats::model.frame(formula, data = data, na.action = stats::na.pass),
    error = function(e) NULL
  )
  if (is.null(mf) || nrow(mf) != nrow(data)) {
    ## A formula the model frame cannot evaluate row-wise: fall back to the columns it names, which
    ## is what the previous implementation always did.
    vars <- intersect(all.vars(formula), names(data))
    return(stats::complete.cases(data[vars]))
  }
  usable <- function(col) {
    if (is.matrix(col)) return(apply(col, 1L, function(r) all(is.finite(r) | !is.numeric(r))))
    if (is.numeric(col)) return(is.finite(col))
    !is.na(col)
  }
  Reduce(`&`, lapply(mf, usable))
}

## Whether the ggplot2 backend is available.
##
## A named function rather than an inline requireNamespace() call, so that the fallback branches are
## reachable in tests on a machine that does have ggplot2: a base function cannot be mocked if the
## package namespace has no binding for it.
has_ggplot2 <- function() requireNamespace("ggplot2", quietly = TRUE)

## Warn when a covariance block cannot support the model it belongs to.
##
## Called from the curve path as well as the tests, so that a design which cannot identify the
## spline is flagged wherever the user enters. Without this, svyrcs_curve() returned a confidence
## band from a rank-deficient covariance without a word.
warn_if_rank_deficient <- function(V, var) {
  k <- ncol(V)
  if (is.null(k) || k < 1L) return(invisible(FALSE))
  r <- tryCatch(qr(V)$rank, error = function(e) k)
  if (r >= k) return(invisible(FALSE))
  warning("the spline covariance for '", var, "' has rank ", r, " but ", k, " columns, so the ",
          "design cannot identify this model. Confidence intervals and tests are not trustworthy; ",
          "use fewer knots, or a design with more primary sampling units.", call. = FALSE)
  invisible(TRUE)
}

## Resolve the denominator degrees of freedom.
##
## The design df is a default, not a law. It is the conservative choice -- it behaves as though every
## covariate varied only between PSUs -- and with few PSUs it matters a great deal. A user analysing
## individual-level covariates, or content with the asymptotic approximation, may reasonably want a
## larger value or Inf. survey::regTermTest() exposes the same choice.
resolve_df <- function(df, degf) {
  if (is.null(df)) {
    if (is.null(degf) || is.na(degf) || (is.finite(degf) && degf <= 0)) {
      stop_svyrcs("could not determine the survey degrees of freedom; pass `df` explicitly, or ",
                  "`df = Inf` for normal-quantile intervals")
    }
    return(degf)
  }
  if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0) {
    stop_svyrcs("`df` must be a single positive number, or Inf for normal-quantile intervals; got ",
                format(df))
  }
  df
}

## Is x a whole number, allowing for the usual floating point slop?
is_count <- function(x, tol = .Machine$double.eps^0.5) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && abs(x - round(x)) < tol
}

## Design-weighted quantiles of a single variable.
##
## survey::svyquantile()'s return shape has changed across versions and also depends on `ci`: with
## ci = FALSE it is a 1 x nq matrix, with ci = TRUE an nq x 4 one, and older survey returned a bare
## numeric. coef() normalises all of them to the estimates alone.
weighted_quantile <- function(var, design, probs) {
  f <- stats::reformulate(var)
  q <- survey::svyquantile(f, design, quantiles = probs, na.rm = TRUE, ci = FALSE)
  out <- as.numeric(stats::coef(q))
  if (length(out) != length(probs) || anyNA(out)) {
    stop_svyrcs("could not compute weighted quantiles of '", var,
                "'; check that it has non-missing values with positive weights")
  }
  stats::setNames(out, paste0(probs * 100, "%"))
}

weighted_mean <- function(var, design) {
  m <- survey::svymean(stats::reformulate(var), design, na.rm = TRUE)
  as.numeric(m)[1L]
}

## Weighted quantities averaged over a list of designs.
##
## With one design these are the plain weighted quantities. With several, they average across
## imputations, which is what fixes the knots and the reference value to a single set of numbers
## before any model is fitted -- without that, the per-imputation coefficient vectors would not be
## estimates of the same parameters and Rubin's rules would not apply.
pooled_weighted_quantile <- function(var, designs, probs) {
  if (length(designs) == 1L) return(weighted_quantile(var, designs[[1L]], probs))
  each <- vapply(designs, function(d) unname(weighted_quantile(var, d, probs)),
                 numeric(length(probs)))
  if (!is.matrix(each)) each <- matrix(each, nrow = length(probs))
  stats::setNames(rowMeans(each), paste0(probs * 100, "%"))
}

pooled_weighted_mean <- function(var, designs) {
  mean(vapply(designs, function(d) weighted_mean(var, d), numeric(1L)))
}

## The exposure values actually used by a fitted model, i.e. after complete-case deletion.
model_variable <- function(fit, var) {
  mf <- stats::model.frame(fit)
  if (!is.null(mf[[var]])) return(mf[[var]])
  ## When the variable only enters through rcs(), the model frame column is the basis matrix and
  ## carries the original values in its first column (which is x itself, unnormalised).
  for (nm in names(mf)) {
    col <- mf[[nm]]
    if (inherits(col, "svyrcs_basis") && identical(attr(col, "var"), var)) return(as.numeric(col[, 1L]))
  }
  NULL
}

## Format a number for printing without scientific notation creeping into reports.
fmt_num <- function(x, digits = 3) {
  formatC(x, digits = digits, format = "g", flag = "#")
}

fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), "NA",
         ifelse(p < 2e-16, "<2e-16", formatC(p, digits = digits, format = "g")))
}
