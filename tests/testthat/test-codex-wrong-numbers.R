## Regressions for the three defects that produced wrong numbers with nothing in the output to show
## it. Each test is built from the reproduction that exposed the bug.

test_that("a transformation that turns finite values non-finite still defines the analytic sample", {
  ## complete.cases() on the raw columns cannot see a NaN that only appears once the formula is
  ## evaluated. Here z <= 0 sits entirely in the upper part of the exposure range, so the rows
  ## log(z) discards are exactly the ones that would drag the outer knot upward.
  set.seed(2)
  n <- 4000
  d <- data.frame(x = runif(n, 10, 60))
  d$z <- ifelse(d$x > 45, -abs(rnorm(n)), abs(rnorm(n)) + 0.5)
  d$y <- rnorm(n)
  d$psu <- rep(1:40, each = n / 40)
  d$strata <- rep(1:20, each = n / 20)
  d$w <- 1
  des <- survey::svydesign(id = ~psu, strata = ~strata, weights = ~w, nest = TRUE, data = d)

  fit <- suppressWarnings(svyrcs(y ~ rcs(x, 4) + log(z), design = des, n = 5))
  used <- rownames(model.frame(fit$model))
  analytic <- des[rownames(d) %in% used, ]

  expect_equal(fit$n, length(used))
  expect_equal(fit$knots,
               signif(unname(coef(survey::svyquantile(~x, analytic, c(0.05, 0.35, 0.65, 0.95),
                                                      na.rm = TRUE, ci = FALSE))), 7))
  expect_equal(fit$ref$value,
               as.numeric(coef(survey::svyquantile(~x, analytic, 0.5, na.rm = TRUE, ci = FALSE))))

  ## and the full-design knots, which the old implementation used, really are different
  full <- signif(unname(coef(survey::svyquantile(~x, des, c(0.05, 0.35, 0.65, 0.95),
                                                 na.rm = TRUE, ci = FALSE))), 7)
  expect_gt(abs(fit$knots[4] - full[4]), 5)
})

test_that("subset= is refused rather than silently applied after knot derivation", {
  des <- nhanes_design()
  expect_error(svyrcs(tchol ~ rcs(bmi, 4), design = des, subset = (age > 40)),
               class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcs(bmi, 4), design = des, subset = (age > 40)),
               "subset\\(design")
})

test_that("a namespaced rcs() freezes its knots for direct model prediction", {
  ## svyrcs() inlines knots itself, so every namespaced test went through a path that never reaches
  ## makepredictcall(). A direct model fit does reach it.
  set.seed(1)
  d <- data.frame(x = runif(300, 1, 100))
  d$y <- 0.02 * d$x + rnorm(300)
  kn <- rcs_knots(d$x, 4)
  nd <- data.frame(x = c(995, 1000, 1005, 1010))
  ref <- rcs_design_matrix(nd$x, kn)

  for (rhs in c("rcs(x, 4)", "svyrcs::rcs(x, 4)", "svyrcs:::rcs(x, 4)")) {
    fit <- lm(stats::as.formula(paste("y ~", rhs)), data = d)
    label <- attr(terms(fit), "term.labels")[1L]
    basis <- model.frame(delete.response(terms(fit)), nd)[[label]]

    expect_equal(attr(basis, "knots"), kn, info = rhs)
    expect_equal(bare(basis), bare(ref), info = rhs)
  }
})

test_that("a pooled imputation fit keeps its designs reachable", {
  des <- mi_design()
  fit <- svyrcs(bmi ~ rcs(age, 4) + tchol, design = des, n = 5)

  dl <- fit_design(fit$model)
  expect_type(dl, "list")
  expect_length(dl, fit$model$m)
  expect_s3_class(dl[[1]], "survey.design")

  ## and the designs are read back from the component fits rather than stored a second time
  expect_null(fit$model$designs)
})

test_that("a named reference on an imputation fit resolves as it did at fitting time", {
  des <- mi_design()
  fit <- svyrcs(bmi ~ rcs(age, 4) + tchol, design = des, n = 5)

  ## the fitting path averages the per-imputation survey-weighted medians
  pooled <- mean(vapply(des$designs, function(d) {
    as.numeric(coef(survey::svyquantile(~age, d, 0.5, na.rm = TRUE, ci = FALSE)))
  }, numeric(1)))
  expect_equal(fit$ref$value, pooled)

  ## asking predict() for the same named reference must land on the same value, so the estimate at
  ## the fitting reference is exactly the null. Before the fix this returned -0.0881, because the
  ## unweighted quantile of the first imputation was used instead.
  expect_equal(predict(fit, x = fit$ref$value, ref = "median")$estimate, 0)
  expect_equal(predict(fit, x = fit$ref$value, ref = "quantile")$estimate, 0)

  ## a different named reference legitimately differs, and a numeric one is unaffected
  expect_true(abs(predict(fit, x = fit$ref$value, ref = "mean")$estimate) > 1e-6)
  expect_equal(predict(fit, x = 50, ref = 50)$estimate, 0)
})

test_that("the pooled fit has not grown", {
  ## Storing the designs alongside the fits pushed object.size() from 47 MB to 59 MB; they are read
  ## back from the component fits instead.
  fit <- svyrcs(bmi ~ rcs(age, 4) + tchol, design = mi_design(), n = 5)
  expect_lt(as.numeric(object.size(fit)) / 1024^2, 50)
})
