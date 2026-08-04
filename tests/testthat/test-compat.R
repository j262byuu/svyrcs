## Regression guard for the ungrouped code path.
##
## Interaction support routes the ungrouped case through the same functions with an identity
## selection matrix and a NULL modifier. These tests pin the behaviour so that a regression in that
## shared path fails loudly instead of quietly shifting numbers.
##
## The numbers were re-baselined once, in 0.5.0, when knots and the reference moved from the design
## as supplied to the analytic sample. `tchol` is missing for 649 rows, so this fixture is exactly
## the case the change is about. Each value below was verified to equal the same quantity computed
## on an independently constructed complete-case design before it was written down; they were not
## copied from the new output.

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
    "  Exposure   bmi, 4 knots at 19.85, 25.27, 29.71, 40.65 (weighted quantiles)",
    "  Sample     9,968 observations, 31 design df",
    "  Reference  bmi = 27.37 (weighted median), Difference = 0",
    "  Dropped    649 rows with missing values",
    "",
    "  Overall association  F = 39.86 on 3 and 31 df,  p = 9.36e-11 ",
    "  Non-linearity        F = 55.21 on 2 and 31 df,  p = 6.07e-11 "
  ))
})

test_that("knots and the reference come from the analytic sample, not the design as supplied", {
  ## The fixture drops 649 rows for missing tchol, so the two differ. Recomputed here from an
  ## independently built complete-case design rather than trusted from the fit.
  fit <- gaussian_fit()
  design <- nhanes_design()
  keep <- stats::complete.cases(nhanes_bmi[c("tchol", "bmi", "age", "sex")])
  analytic <- design[keep, ]

  kn <- as.numeric(coef(survey::svyquantile(~bmi, analytic, c(0.05, 0.35, 0.65, 0.95),
                                            na.rm = TRUE, ci = FALSE)))
  expect_equal(fit$knots, signif(kn, 7))
  expect_equal(fit$ref$value,
               as.numeric(coef(survey::svyquantile(~bmi, analytic, 0.5, na.rm = TRUE, ci = FALSE))))
  expect_equal(range(fit$curve$x),
               as.numeric(coef(survey::svyquantile(~bmi, analytic, c(0.01, 0.99),
                                                   na.rm = TRUE, ci = FALSE))))

  ## and they really are different from the full-design values, so the test has teeth
  kn_full <- as.numeric(coef(survey::svyquantile(~bmi, design, c(0.05, 0.35, 0.65, 0.95),
                                                 na.rm = TRUE, ci = FALSE)))
  expect_false(isTRUE(all.equal(signif(kn, 7), signif(kn_full, 7))))
})

test_that("ungrouped estimates are unchanged", {
  fit <- gaussian_fit()
  p <- predict(fit, x = c(18, 25, 30, 40))
  expect_equal(round(p$estimate, 6), c(-20.466448, -2.527119, 0.669398, -3.563271))
  expect_equal(round(p$se, 6), c(2.079350, 0.414612, 0.595546, 1.188520))
  expect_equal(round(fit$knots, 3), c(19.85, 25.27, 29.71, 40.65))
  expect_equal(round(fit$ref$value, 3), 27.37)
  expect_equal(round(fit$tests$overall$F, 4), 39.8636)
  expect_equal(round(fit$tests$nonlinear$F, 4), 55.2097)
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
