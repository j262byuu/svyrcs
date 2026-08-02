## Regression guard for the ungrouped code path.
##
## Interaction support routes the ungrouped case through the same functions with an identity
## selection matrix and a NULL modifier. These tests pin the 0.1.0 behaviour so that a regression in
## that shared path fails loudly instead of quietly shifting numbers.

test_that("the ungrouped fit has no group column and the 0.1.0 object shape", {
  fit <- gaussian_fit()
  expect_null(fit$groups)
  expect_false("group" %in% names(fit$curve))
  expect_named(fit$tests, c("overall", "nonlinear"))
  expect_named(fit$ref, c("value", "method"))
  expect_equal(nrow(fit$curve), 200L)
  expect_named(as.data.frame(fit$curve),
               c("x", "estimate", "conf.low", "conf.high", "se"))
})

test_that("ungrouped print output is unchanged", {
  fit <- gaussian_fit()
  out <- capture.output(print(fit))
  expect_equal(out, c(
    "Restricted cubic spline on a complex survey design",
    "",
    "  Model      survey-weighted GLM (gaussian, identity link)",
    "  Outcome    tchol",
    "  Exposure   bmi, 4 knots at 19.78, 25.22, 29.71, 40.67 (weighted quantiles)",
    "  Sample     9,968 observations, 31 design df",
    "  Reference  bmi = 27.35 (weighted median), Difference = 0",
    "  Dropped    649 rows with missing values",
    "",
    "  Overall association  F = 39.87 on 3 and 31 df,  p = 9.34e-11 ",
    "  Non-linearity        F = 55.21 on 2 and 31 df,  p = 6.07e-11 "
  ))
})

test_that("ungrouped estimates are unchanged", {
  fit <- gaussian_fit()
  p <- predict(fit, x = c(18, 25, 30, 40))
  expect_equal(round(p$estimate, 6), c(-20.458087, -2.505363, 0.692455, -3.530227))
  expect_equal(round(p$se, 6), c(2.081394, 0.412191, 0.597615, 1.188749))
  expect_equal(round(fit$knots, 3), c(19.78, 25.22, 29.71, 40.67))
  expect_equal(round(fit$ref$value, 3), 27.35)
  expect_equal(round(fit$tests$overall$F, 4), 39.8715)
  expect_equal(round(fit$tests$nonlinear$F, 4), 55.2066)
})

test_that("the ungrouped Cox fit is unchanged", {
  skip_if_not_installed("survival")
  fit <- cox_fit()
  expect_equal(fit$n, 10617L)
  expect_equal(fit$nevents, 1893L)
  expect_equal(fit$degf, 31)
  expect_equal(round(fit$tests$overall$F, 4), 48.0731)
  expect_equal(round(predict(fit, x = c(18, 40))$estimate, 6), c(2.785167, 1.499309))
})
