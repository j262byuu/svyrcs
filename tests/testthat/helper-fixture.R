## Shared fixtures, built once per test run.
##
## Everything is fitted on the real shipped NHANES design rather than on simulated data, so the tests
## exercise the same code paths users will: a stratified multistage design with 31 degrees of
## freedom, missing covariates, and a genuinely non-linear exposure.

nhanes_design <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- survey::svydesign(
        id = ~psu, strata = ~strata, weights = ~weight,
        nest = TRUE, data = svyrcs::nhanes_bmi
      )
    }
    cached
  }
})

cox_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- svyrcs(survival::Surv(time, event) ~ rcs(bmi, 4) + age + sex,
                        design = nhanes_design())
    }
    cached
  }
})

gaussian_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- svyrcs(tchol ~ rcs(bmi, 4) + age + sex, design = nhanes_design())
    }
    cached
  }
})

bare_glm <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- survey::svyglm(tchol ~ rcs(bmi, 4) + age + sex, design = nhanes_design())
    }
    cached
  }
})

## Strip a basis down to a plain numeric matrix, so comparisons are about numbers rather than
## the knot / class attributes that rcs() deliberately attaches.
bare <- function(b) {
  b <- unclass(b)
  attributes(b) <- list(dim = dim(b))
  b
}
