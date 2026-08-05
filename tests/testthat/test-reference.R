test_that("the median reference is the survey-weighted median, not the sample median", {
  fit <- gaussian_fit()
  ## the analytic sample, since that is what the fit is based on
  design <- nhanes_design()[stats::complete.cases(nhanes_bmi[c("tchol", "bmi", "age", "sex")]), ]
  wq <- as.numeric(coef(survey::svyquantile(~bmi, design, 0.5, na.rm = TRUE, ci = FALSE)))
  expect_equal(fit$ref$value, wq)
  expect_equal(fit$ref$method, "weighted median")
  ## and it genuinely differs from the unweighted one
  expect_false(isTRUE(all.equal(wq, median(nhanes_bmi$bmi))))
})

test_that("the mean reference is the survey-weighted mean", {
  fit <- svyrcs(tchol ~ rcspline(bmi, 4) + age, design = nhanes_design(), ref = "mean")
  analytic <- nhanes_design()[stats::complete.cases(nhanes_bmi[c("tchol", "bmi", "age")]), ]
  wm <- as.numeric(survey::svymean(~bmi, analytic, na.rm = TRUE))[1L]
  expect_equal(fit$ref$value, wm)
  expect_equal(fit$ref$method, "weighted mean")
})

test_that("an arbitrary weighted quantile can be used as the reference", {
  fit <- svyrcs(tchol ~ rcspline(bmi, 4) + age, design = nhanes_design(),
                ref = "quantile", ref_prob = 0.25)
  analytic <- nhanes_design()[stats::complete.cases(nhanes_bmi[c("tchol", "bmi", "age")]), ]
  wq <- as.numeric(coef(survey::svyquantile(~bmi, analytic, 0.25, na.rm = TRUE, ci = FALSE)))
  expect_equal(fit$ref$value, wq)
  expect_match(fit$ref$method, "25th percentile")
})

test_that("a numeric reference is used as given", {
  fit <- svyrcs(tchol ~ rcspline(bmi, 4) + age, design = nhanes_design(), ref = 25)
  expect_equal(fit$ref$value, 25)
  expect_equal(fit$ref$method, "user-specified")
  expect_equal(predict(fit, x = 25)$estimate, 0)
})

test_that("min and max references land on the extremes of the curve", {
  skip_if_not_installed("survival")
  fmin <- svyrcs(survival::Surv(time, event) ~ rcspline(bmi, 4) + age + sex,
                 design = nhanes_design(), ref = "min", n = 500)
  fmax <- svyrcs(survival::Surv(time, event) ~ rcspline(bmi, 4) + age + sex,
                 design = nhanes_design(), ref = "max", n = 500)

  expect_match(fmin$ref$method, "^minimum-risk point")
  expect_match(fmax$ref$method, "^maximum-risk point")
  ## the maximum-risk BMI here is the low end of the plotting range, and the label says so
  expect_match(fmax$ref$method, "at the range boundary")
  expect_false(grepl("range boundary", fmin$ref$method))

  ## the curve anchored at its own minimum cannot dip below the null
  expect_gte(min(fmin$curve$estimate), 1 - 1e-8)
  expect_lte(max(fmax$curve$estimate), 1 + 1e-8)

  ## the located points agree with the argmin/argmax of an independently referenced curve
  ref_curve <- predict(fmin, x = seq(min(fmin$curve$x), max(fmin$curve$x), length.out = 500),
                       ref = 25)
  expect_equal(fmin$ref$value, ref_curve$x[which.min(ref_curve$estimate)], tolerance = 0.05)
  expect_equal(fmax$ref$value, ref_curve$x[which.max(ref_curve$estimate)], tolerance = 0.05)
})

test_that("the reference choice shifts the curve without changing its shape", {
  f1 <- gaussian_fit()
  f2 <- svyrcs(tchol ~ rcspline(bmi, 4) + age + sex, design = nhanes_design(), ref = 30)
  ## for a difference measure, changing the reference is a pure vertical shift
  shift <- f1$curve$estimate - f2$curve$estimate
  expect_equal(diff(range(shift)), 0, tolerance = 1e-8)
  ## and the tests are unaffected, since they are about the coefficients
  expect_equal(f1$tests$overall$F, f2$tests$overall$F)
})

test_that("invalid references are rejected", {
  design <- nhanes_design()
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = design, ref = 999),
               class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = design, ref = "middle"),
               class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = design, ref = c(20, 30)),
               class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = design, ref = "quantile", ref_prob = 1.2),
               class = "svyrcs_error")
})
