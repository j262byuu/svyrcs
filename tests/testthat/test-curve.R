test_that("the estimate at the reference is exactly null with zero standard error", {
  skip_if_not_installed("survival")
  fit <- cox_fit()
  at_ref <- predict(fit, x = fit$ref$value)
  expect_identical(at_ref$estimate, 1)
  expect_identical(at_ref$se, 0)
  expect_identical(at_ref$conf.low, 1)
  expect_identical(at_ref$conf.high, 1)

  g <- gaussian_fit()
  at_ref_g <- predict(g, x = g$ref$value)
  expect_identical(at_ref_g$estimate, 0)
  expect_identical(at_ref_g$se, 0)
})

test_that("the contrast reproduces predict(type = 'lp') differences", {
  skip_if_not_installed("survival")
  fit <- cox_fit()
  design <- nhanes_design()

  grid <- fit$curve$x
  nd <- design$variables[rep(1L, length(grid)), ]
  nd$bmi <- grid
  lp <- as.numeric(predict(fit$model, newdata = nd, type = "lp"))

  nd0 <- nd[1L, ]
  nd0$bmi <- fit$ref$value
  lp0 <- as.numeric(predict(fit$model, newdata = nd0, type = "lp"))

  expect_equal(log(fit$curve$estimate), lp - lp0, tolerance = 1e-8)
})

test_that("covariate values do not affect the curve, because they cancel", {
  skip_if_not_installed("survival")
  fit <- cox_fit()
  design <- nhanes_design()

  ## Two different covariate profiles must give the same contrast.
  contrast_from <- function(row_index) {
    nd <- design$variables[rep(row_index, 2L), ]
    nd$bmi <- c(22, 34)
    lp <- as.numeric(predict(fit$model, newdata = nd, type = "lp"))
    lp[2L] - lp[1L]
  }
  young <- which(design$variables$age < 30 & design$variables$sex == "Male")[1L]
  old <- which(design$variables$age > 70 & design$variables$sex == "Female")[1L]
  expect_equal(contrast_from(young), contrast_from(old), tolerance = 1e-10)

  p <- predict(fit, x = c(22, 34), ref = 22)
  expect_equal(log(p$estimate[2L]), contrast_from(young), tolerance = 1e-8)
})

test_that("a covariate whose name contains the exposure name is not absorbed", {
  design <- nhanes_design()
  design$variables$bmi_over30 <- as.numeric(design$variables$bmi > 30)
  fit <- svyrcs(tchol ~ rcs(bmi, 4) + bmi_over30 + age, design = design)

  term <- find_rcs_term(fit$model)
  idx <- spline_coef_index(fit$model, term)
  expect_length(idx, fit$nk - 1L)
  expect_false(any(grepl("bmi_over30", names(coef(fit$model))[idx], fixed = TRUE)))
})

test_that("svyrcs_curve works on a model the user fitted themselves", {
  m <- bare_glm()
  cv <- svyrcs_curve(m, "bmi", n = 25)
  expect_s3_class(cv, "svyrcs_curve")
  expect_equal(nrow(cv), 25L)
  expect_equal(attr(cv, "measure"), "Difference")
  expect_equal(attr(cv, "degf"), survey::degf(nhanes_design()))
  expect_true(all(is.finite(cv$estimate)))
  expect_true(all(cv$conf.low <= cv$estimate))
  expect_true(all(cv$conf.high >= cv$estimate))
})

test_that("`at` overrides the grid and `range` bounds it", {
  m <- bare_glm()
  at <- c(20, 25, 30)
  cv <- svyrcs_curve(m, "bmi", at = at)
  expect_equal(cv$x, at)

  cv2 <- svyrcs_curve(m, "bmi", range = c(20, 30), n = 11)
  expect_equal(range(cv2$x), c(20, 30))
  expect_equal(nrow(cv2), 11L)
})

test_that("the confidence level uses a t quantile on the design degrees of freedom", {
  m <- bare_glm()
  cv <- svyrcs_curve(m, "bmi", at = c(20, 35), level = 0.95)
  degf <- survey::degf(nhanes_design())
  crit <- qt(0.975, degf)
  expect_equal(cv$conf.high - cv$estimate, crit * cv$se, tolerance = 1e-10)

  wide <- svyrcs_curve(m, "bmi", at = c(20, 35), level = 0.99)
  expect_true(all(wide$conf.high - wide$estimate > cv$conf.high - cv$estimate))
  ## a normal quantile would be visibly narrower with only 31 df
  expect_gt(crit, qnorm(0.975))
})

test_that("curve estimation rejects bad input", {
  m <- bare_glm()
  expect_error(svyrcs_curve(m, "nosuchvar"), class = "svyrcs_error")
  expect_error(svyrcs_curve(m, "bmi", at = character(0)), class = "svyrcs_error")
  expect_error(svyrcs_curve(m, "bmi", n = 1), class = "svyrcs_error")
  expect_error(svyrcs_curve(m, "bmi", level = 1.5), class = "svyrcs_error")
  expect_error(svyrcs_curve(lm(tchol ~ bmi, data = nhanes_bmi), "bmi"), class = "svyrcs_error")
})
