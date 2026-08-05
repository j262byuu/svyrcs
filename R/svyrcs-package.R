#' @keywords internal
#' @aliases svyrcs-package
#'
#' @section Getting started:
#' Build a survey design with [survey::svydesign()], then call [svyrcs()] with a model formula
#' containing one [rcspline()] term:
#'
#' ```r
#' design <- svydesign(id = ~psu, strata = ~strata, weights = ~weight,
#'                     nest = TRUE, data = nhanes_bmi)
#' fit <- svyrcs(Surv(time, event) ~ rcspline(bmi, 4) + age + sex, design = design)
#' summary(fit)
#' plot(fit)
#' ```
#'
#' @section Using this package alongside `rms`:
#' The spline function here is [rcspline()]. It was called `rcs()` before 0.7.0, which masked
#' `rms::rcs()` whenever both packages were attached; the name was changed because that masking could
#' not be made safe. A `cph()`, `lrm()` or `ols()` model written with a bare `rcs()` would silently
#' pick up this basis, fit correctly, and then lose the `Nonlinear` row from `anova()` -- with no
#' error and no warning, and that row is usually the reason for fitting a spline at all.
#'
#' Since the rename there is no collision, and no attach order to remember. Use `rms::rcs()` in `rms`
#' models and [rcspline()] here. Note that the converse is not symmetric: passing [rcspline()] into an
#' `rms` model still fits and still loses the `Nonlinear` test, because the basis carries this
#' package's knot attributes rather than `rms`'s design attributes.
"_PACKAGE"

## Every generic this package registers an S3 method for must be imported here, or the namespace
## cannot be loaded with only base attached.
#' @importFrom stats coef vcov nobs formula predict makepredictcall
#' @importFrom stats qt pf pchisq quantile terms model.frame
#' @importFrom stats setNames reformulate gaussian median family
#' @importFrom survey svydesign svyglm svycoxph degf svyquantile svymean
#' @importFrom rlang .data
NULL
