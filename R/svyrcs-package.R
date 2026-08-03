#' @keywords internal
#' @aliases svyrcs-package
#'
#' @section Getting started:
#' Build a survey design with [survey::svydesign()], then call [svyrcs()] with a model formula
#' containing one [rcs()] term:
#'
#' ```r
#' design <- svydesign(id = ~psu, strata = ~strata, weights = ~weight,
#'                     nest = TRUE, data = nhanes_bmi)
#' fit <- svyrcs(Surv(time, event) ~ rcs(bmi, 4) + age + sex, design = design)
#' summary(fit)
#' plot(fit)
#' ```
#'
#' @section Note on masking:
#' `svyrcs` exports a function called `rcs()`, which masks `rms::rcs()` when both packages are
#' attached. The two produce numerically identical bases (agreement to about 1e-14), so the masking
#' does not change any result. It does add knot storage, which is what makes [predict()] on new data
#' safe. See [rcs()].
"_PACKAGE"

## Every generic this package registers an S3 method for must be imported here, or the namespace
## cannot be loaded with only base attached.
#' @importFrom stats coef vcov nobs formula predict makepredictcall
#' @importFrom stats qt pf pchisq quantile terms model.frame
#' @importFrom stats setNames reformulate gaussian median family
#' @importFrom survey svydesign svyglm svycoxph degf svyquantile svymean
#' @importFrom rlang .data
NULL
