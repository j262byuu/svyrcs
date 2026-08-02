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
