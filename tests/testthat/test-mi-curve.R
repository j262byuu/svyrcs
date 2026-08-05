mi_gauss <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- svyrcs(bmi ~ rcspline(age, 4) + sex + tchol + hdl, design = mi_design())
    }
    cached
  }
})

test_that("an imputed fit carries per-point degrees of freedom and fmi", {
  fit <- mi_gauss()
  expect_true(all(c("df", "fmi") %in% names(fit$curve)))
  expect_equal(nrow(fit$curve), 200L)
  expect_true(all(is.finite(fit$curve$df)))
  expect_true(all(fit$curve$fmi >= 0 & fit$curve$fmi < 1))
  expect_equal(fit$imputations$m, 5L)
})

test_that("the degrees of freedom never exceed the design degrees of freedom", {
  fit <- mi_gauss()
  expect_true(all(fit$curve$df <= fit$degf + 1e-9))
  expect_gt(max(fit$curve$fmi), 0)
  ## and the points with more missing information have fewer degrees of freedom
  expect_lt(cor(fit$curve$fmi, fit$curve$df), 0)
})

test_that("intervals are wider than the design degrees of freedom alone would give", {
  fit <- mi_gauss()
  cv <- fit$curve
  naive_half <- qt(0.975, fit$degf) * cv$se
  mi_half <- cv$conf.high - cv$estimate
  ## equal only where fmi is 0; wider everywhere else
  informative <- cv$fmi > 1e-8
  expect_true(any(informative))
  expect_true(all(mi_half[informative] > naive_half[informative]))
  expect_equal(mi_half[!informative], naive_half[!informative])
})

test_that("our degrees of freedom are smaller than mitools' uncorrected rule", {
  skip_if_not_installed("mitools")
  fit <- mi_gauss()
  m <- fit$model$m
  term <- find_rcs_term(fit$model)
  L <- selection_matrix(names(coef(fit$model)), spline_colnames(term$label, term$nk), NULL)

  x <- fit$curve$x
  dB <- sweep(rcs_design_matrix(x, term$knots), 2,
              as.numeric(rcs_design_matrix(fit$ref$value, term$knots)))
  M <- dB %*% L
  u <- rowSums((M %*% fit$model$Ubar) * M)
  b <- rowSums((M %*% fit$model$B) * M)

  informative <- b > 0
  r <- (1 + 1 / m) * b[informative] / u[informative]
  mitools_df <- (m - 1) * (1 + 1 / r)^2
  expect_true(all(fit$curve$df[informative] < mitools_df))
  ## and mitools' rule really is the unusable one here
  expect_gt(max(mitools_df), 10 * fit$degf)
})

test_that("the joint tests use an adjusted denominator, bounded by the design df", {
  fit <- mi_gauss()
  expect_lte(fit$tests$overall$df2, fit$degf + 1e-9)
  expect_lte(fit$tests$nonlinear$df2, fit$degf + 1e-9)
  expect_true(is.finite(fit$tests$overall$p_F))
  expect_gte(fit$tests$overall$fmi, 0)
  ## the statistic is D1's, on (1 + r) Ubar rather than the pooled total covariance
  term <- find_rcs_term(fit$model)
  cols <- spline_colnames(term$label, term$nk)
  b <- coef(fit$model)[cols]
  d1 <- d1_test(b, fit$model$Ubar[cols, cols], fit$model$B[cols, cols], fit$model$m, fit$degf)
  expect_equal(fit$tests$overall$F, d1$F, tolerance = 1e-10)
  expect_equal(fit$tests$overall$df2, d1$v, tolerance = 1e-10)
})

test_that("with nothing imputed the degrees of freedom take the formula's zero-information value", {
  ## The Barnard-Rubin and Reiter formulas do not return the complete-data df when the
  ## between-imputation variance is zero; they return (nu + 1) nu / (nu + 3), about 6% lower at
  ## nu = 31. Following them literally is the point of using a named rule.
  ## m = 3 and a two-parameter non-linear block give k(m - 1) = 4, where Reiter's expansion is
  ## undefined and the fallback warns. That is asserted separately below.
  fit <- suppressWarnings(svyrcs(bmi ~ rcspline(age, 4) + sex, design = mi_design_degenerate()))
  vstar <- (fit$degf + 1) / (fit$degf + 3) * fit$degf
  expect_equal(unique(fit$curve$df), vstar)
  expect_identical(unique(fit$curve$fmi), 0)
  expect_equal(fit$tests$overall$df2, vstar)
  expect_equal(fit$tests$overall$fmi, 0)
  expect_lt(vstar, fit$degf)
})

test_that("a block too small for the small-sample correction warns and falls back", {
  ## Reiter's expansion contains 1/(k(m - 1) - 4). With m = 3 and a 2-parameter block that is a
  ## division by zero; mitml returns 4 or NaN there, which is not a usable reference distribution.
  w <- capture_warnings(fit <- svyrcs(bmi ~ rcspline(age, 4) + sex, design = mi_design_degenerate()))
  expect_true(any(grepl("k\\(m - 1\\) > 4", w)))
  expect_true(any(grepl("increase the number of imputations", w)))
  expect_equal(fit$tests$nonlinear$df2, fit$degf)
  expect_true(is.finite(fit$tests$nonlinear$p_F))

  ## with more imputations the correction applies and no warning is raised
  big <- svyrcs(bmi ~ rcspline(age, 4) + sex,
                design = survey::svydesign(id = ~psu, strata = ~strata, weights = ~weight,
                                           nest = TRUE,
                                           data = mitools::imputationList(rep(list(nhanes_bmi), 6))))
  expect_lt(big$tests$nonlinear$df2, big$degf)
})

test_that("an ordinary fit gains no df or fmi columns", {
  fit <- gaussian_fit()
  expect_false("df" %in% names(fit$curve))
  expect_false("fmi" %in% names(fit$curve))
  expect_null(fit$imputations)
  expect_null(fit$tests$overall$fmi)
})
