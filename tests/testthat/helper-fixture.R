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
      cached <<- svyrcs(survival::Surv(time, event) ~ rcspline(bmi, 4) + age + sex,
                        design = nhanes_design())
    }
    cached
  }
})

gaussian_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- svyrcs(tchol ~ rcspline(bmi, 4) + age + sex, design = nhanes_design())
    }
    cached
  }
})

bare_glm <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- survey::svyglm(tchol ~ rcspline(bmi, 4) + age + sex, design = nhanes_design())
    }
    cached
  }
})

## Fit a model with an explicit-knot rcspline() term directly, bypassing svyrcs(), so the curve and
## modifier machinery can be tested independently of the formula surgery in svyrcs().
##
## Knots must be inlined: survey::svyglm() and svycoxph() reject a formula referring to any symbol
## that is not a column of the design.
int_model <- function(rhs, response = "tchol", cox = FALSE, data_extra = NULL, nk = 4) {
  design <- nhanes_design()
  if (!is.null(data_extra)) {
    for (nm in names(data_extra)) design$variables[[nm]] <- data_extra[[nm]]
  }
  kn <- rcs_knots(design$variables$bmi, nk)
  txt <- paste(response, "~", rhs)
  f <- stats::as.formula(do.call(substitute, list(str2lang(txt), list(KN = kn))))
  if (cox) survey::svycoxph(f, design = design) else survey::svyglm(f, design = design)
}

## Strip a basis down to a plain numeric matrix, so comparisons are about numbers rather than
## the knot / class attributes that rcspline() deliberately attaches.
bare <- function(b) {
  b <- unclass(b)
  attributes(b) <- list(dim = dim(b))
  b
}

## A replicate design whose vcov() comes back with NaN entries.
##
## Negative `rscales` are not a fabrication: successive-difference replication and some Fay variants
## legitimately carry negative combining coefficients, and survey::svrepdesign() accepts them. The
## variance routine takes sqrt(rscales), so the covariance picks up NaN while the coefficients stay
## finite -- which is exactly the shape that used to slip past the finiteness checks.
nan_vcov_design <- function(frac = 0.4, seed = 4) {
  jk <- survey::as.svrepdesign(nhanes_design(), type = "JKn")
  rs <- jk$rscales
  set.seed(seed)
  rs[sample(seq_along(rs), floor(frac * length(rs)))] <-
    -rs[sample(seq_along(rs), floor(frac * length(rs)))]
  survey::svrepdesign(data = jk$variables, repweights = as.matrix(jk$repweights),
                      weights = stats::weights(jk, "sampling"), type = "other",
                      scale = jk$scale, rscales = rs, mse = TRUE, combined.weights = FALSE)
}

## A genuine replicate-weight design. The suite exercised these for the first time in 0.6.3, and the
## Cox path on one was broken until 0.6.4.
rep_design <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) cached <<- survey::as.svrepdesign(nhanes_design(), type = "JKn")
    cached
  }
})
