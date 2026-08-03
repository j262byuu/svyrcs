test_that("the overall F test reproduces survey::regTermTest", {
  skip_if_not_installed("survival")
  fit <- cox_fit()
  rt <- survey::regTermTest(fit$model, reformulate(fit$term_label), df = fit$degf)
  expect_equal(as.numeric(rt$Ftest), fit$tests$overall$F, tolerance = 1e-8)
  expect_equal(as.numeric(rt$df), fit$tests$overall$df1)
  expect_equal(as.numeric(rt$p), fit$tests$overall$p_F, tolerance = 1e-8)
})

test_that("the overall F test reproduces regTermTest for a glm too", {
  fit <- gaussian_fit()
  rt <- survey::regTermTest(fit$model, reformulate(fit$term_label), df = fit$degf)
  expect_equal(as.numeric(rt$Ftest), fit$tests$overall$F, tolerance = 1e-8)
})

test_that("F and chi-square are consistent, and F is the more conservative", {
  fit <- gaussian_fit()
  for (tt in list(fit$tests$overall, fit$tests$nonlinear)) {
    expect_equal(tt$F, tt$chisq / tt$df1)
    expect_equal(tt$p_chisq, pchisq(tt$chisq, tt$df1, lower.tail = FALSE))
    expect_equal(tt$p_F, pf(tt$F, tt$df1, tt$df2, lower.tail = FALSE))
    ## with finite design df the F reference distribution has heavier tails
    expect_gt(tt$p_F, tt$p_chisq)
  }
})

test_that("the non-linearity test drops exactly the linear term", {
  fit <- gaussian_fit()
  expect_equal(fit$tests$overall$df1, fit$nk - 1L)
  expect_equal(fit$tests$nonlinear$df1, fit$nk - 2L)
  expect_equal(fit$tests$overall$df2, fit$degf)
  expect_equal(fit$tests$nonlinear$df2, fit$degf)
})

test_that("the non-linearity test matches a manual Wald test on the non-linear columns", {
  fit <- gaussian_fit()
  term <- find_rcs_term(fit$model)
  idx <- spline_coef_index(fit$model, term)[-1L]
  b <- coef(fit$model)[idx]
  V <- vcov(fit$model)[idx, idx, drop = FALSE]
  chisq <- as.numeric(t(b) %*% solve(V) %*% b)
  expect_equal(fit$tests$nonlinear$chisq, chisq, tolerance = 1e-8)
})

test_that("3 knots leave a single non-linear column", {
  fit <- svyrcs(tchol ~ rcs(bmi, 3) + age, design = nhanes_design())
  expect_equal(fit$tests$overall$df1, 2L)
  expect_equal(fit$tests$nonlinear$df1, 1L)
  expect_true(is.finite(fit$tests$nonlinear$p_F))
})

test_that("a singular spline block warns and reports the test as unavailable", {
  b <- c(1, 2)
  V <- matrix(c(1, 1, 1, 1), 2, 2)  # exactly singular
  expect_warning(out <- wald_block(b, V, 30, "overall association"), "singular")
  ## warn rather than stop, so that svyrcs() and svyrcs_curve() behave the same way on the same data
  expect_true(is.na(out$chisq))
  expect_true(is.na(out$p_F))
  expect_equal(out$df1, 2L)
})
