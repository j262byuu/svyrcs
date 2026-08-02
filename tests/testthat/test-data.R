test_that("nhanes_bmi matches its documented contract", {
  expect_s3_class(nhanes_bmi, "data.frame")
  expect_equal(nrow(nhanes_bmi), 10617L)
  expect_named(nhanes_bmi,
               c("seqn", "cycle", "psu", "strata", "weight", "age", "sex", "race", "bmi",
                 "tchol", "hdl", "glucose", "statin_use", "time", "event", "high_chol"))
  expect_equal(levels(nhanes_bmi$cycle), c("2005-2006", "2007-2008"))
  expect_equal(levels(nhanes_bmi$sex), c("Male", "Female"))
  expect_equal(nlevels(nhanes_bmi$race), 4L)
})

test_that("the survey design structure is intact", {
  expect_false(anyNA(nhanes_bmi[c("psu", "strata", "weight", "bmi", "age", "sex", "race",
                                  "time", "event")]))
  expect_true(all(nhanes_bmi$weight > 0))
  expect_true(all(nhanes_bmi$time > 0))
  expect_true(all(nhanes_bmi$event %in% c(0L, 1L)))

  design <- nhanes_design()
  expect_equal(survey::degf(design), 31)
  ## every stratum must contain more than one PSU, or the variance is not estimable
  psu_per_stratum <- tapply(nhanes_bmi$psu, nhanes_bmi$strata, function(x) length(unique(x)))
  expect_true(all(psu_per_stratum > 1L))
})

test_that("the BMI-mortality curve reproduces the expected J shape", {
  skip_if_not_installed("survival")
  fit <- cox_fit()

  ## a real, strongly non-linear association
  expect_lt(fit$tests$overall$p_F, 1e-6)
  expect_lt(fit$tests$nonlinear$p_F, 1e-6)

  ## elevated risk at low BMI, a minimum in the mid-to-high 20s, rising again at high BMI
  p <- predict(fit, x = c(18, 27, 40))
  expect_gt(p$estimate[1L], 1.5)
  expect_lt(p$estimate[2L], p$estimate[1L])
  expect_gt(p$estimate[3L], p$estimate[2L])

  low <- fit$curve$x[which.min(fit$curve$estimate)]
  expect_gt(low, 24)
  expect_lt(low, 32)
})
