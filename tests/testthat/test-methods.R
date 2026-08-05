test_that("print reports the model, the exposure, the sample and both tests", {
  skip_if_not_installed("survival")
  out <- capture.output(print(cox_fit()))
  expect_true(any(grepl("Cox proportional hazards", out)))
  expect_true(any(grepl("bmi, 4 knots", out)))
  expect_true(any(grepl("weighted quantiles", out)))
  expect_true(any(grepl("10,617 observations", out)))
  expect_true(any(grepl("1,893 events", out)))
  expect_true(any(grepl("31 design df", out)))
  expect_true(any(grepl("Overall association", out)))
  expect_true(any(grepl("Non-linearity", out)))
})

test_that("print returns its argument invisibly", {
  fit <- gaussian_fit()
  expect_invisible(print(fit))
  expect_identical(withVisible(print(fit))$value, fit)
})

test_that("summary adds the curve shape and the significant ranges", {
  fit <- gaussian_fit()
  s <- summary(fit)
  expect_s3_class(s, "summary.svyrcs")
  expect_named(s$features, c("min", "max", "significant"))
  expect_equal(s$features$min$estimate, min(fit$curve$estimate))
  expect_equal(s$features$max$estimate, max(fit$curve$estimate))

  out <- capture.output(print(s))
  expect_true(any(grepl("Curve shape", out)))
  expect_true(any(grepl("Grid lowest", out)))
  expect_true(any(grepl("Grid highest", out)))
})

test_that("the significant ranges really are where the band excludes the null", {
  fit <- gaussian_fit()
  sig <- summary(fit)$features$significant
  skip_if(nrow(sig) == 0)
  cv <- fit$curve
  for (i in seq_len(nrow(sig))) {
    inside <- cv$x >= sig$from[i] & cv$x <= sig$to[i]
    if (sig$direction[i] == "above") {
      expect_true(all(cv$conf.low[inside] > fit$null))
    } else {
      expect_true(all(cv$conf.high[inside] < fit$null))
    }
  }
  ## and nothing significant was missed
  flagged <- rep(FALSE, nrow(cv))
  for (i in seq_len(nrow(sig))) flagged <- flagged | (cv$x >= sig$from[i] & cv$x <= sig$to[i])
  excl <- cv$conf.low > fit$null | cv$conf.high < fit$null
  expect_equal(flagged, excl)
})

test_that("summary reports no significant range when the band always covers the null", {
  ## an exposure with no association: a curve whose band straddles the null everywhere
  design <- nhanes_design()
  set.seed(11)
  design$variables$noise_outcome <- rnorm(nrow(design$variables))
  fit <- svyrcs(noise_outcome ~ rcspline(bmi, 4) + age, design = design)
  s <- summary(fit)
  if (nrow(s$features$significant) == 0L) {
    expect_output(print(s), "includes")
  } else {
    succeed()
  }
})

test_that("predict evaluates at chosen values and matches the curve", {
  fit <- gaussian_fit()
  at <- fit$curve$x[c(1L, 40L, 120L, 200L)]
  p <- predict(fit, x = at)
  expect_equal(p$x, at)
  expect_equal(p$estimate, fit$curve$estimate[c(1L, 40L, 120L, 200L)])
  expect_equal(p$conf.low, fit$curve$conf.low[c(1L, 40L, 120L, 200L)])
})

test_that("predict can re-reference without refitting", {
  fit <- gaussian_fit()
  p_ref25 <- predict(fit, x = c(20, 30), ref = 25)
  expect_equal(predict(fit, x = 25, ref = 25)$estimate, 0)
  ## for a difference measure, re-referencing is subtraction
  p_orig <- predict(fit, x = c(20, 30, 25))
  expect_equal(p_ref25$estimate, p_orig$estimate[1:2] - p_orig$estimate[3], tolerance = 1e-10)
})

test_that("predict honours a different confidence level", {
  fit <- gaussian_fit()
  p95 <- predict(fit, x = 35)
  p99 <- predict(fit, x = 35, level = 0.99)
  expect_equal(p95$se, p99$se)
  expect_gt(p99$conf.high, p95$conf.high)
})

test_that("predict rejects bad input", {
  fit <- gaussian_fit()
  expect_error(predict(fit), class = "svyrcs_error")
  expect_error(predict(fit, x = "30"), class = "svyrcs_error")
  expect_error(predict(fit, x = numeric(0)), class = "svyrcs_error")
})

test_that("the generic accessors delegate to the fitted model", {
  fit <- gaussian_fit()
  expect_identical(coef(fit), coef(fit$model))
  expect_identical(vcov(fit), vcov(fit$model))
  expect_identical(nobs(fit), fit$n)
  expect_identical(formula(fit), fit$formula)
})

test_that("as.data.frame strips the curve's attributes", {
  fit <- gaussian_fit()
  d <- as.data.frame(fit$curve)
  expect_s3_class(d, "data.frame")
  expect_false(inherits(d, "svyrcs_curve"))
  expect_null(attr(d, "measure"))
  expect_named(d, c("x", "estimate", "conf.low", "conf.high", "se"))
})
