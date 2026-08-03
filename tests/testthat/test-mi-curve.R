mi_gauss <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- svyrcs(bmi ~ rcs(age, 4) + sex + tchol + hdl, design = mi_design())
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
  ## the statistic itself is a plain Wald on the pooled covariance
  term <- find_rcs_term(fit$model)
  cols <- spline_colnames(term$label, term$nk)
  b <- coef(fit$model)[cols]
  V <- vcov(fit$model)[cols, cols, drop = FALSE]
  expect_equal(fit$tests$overall$chisq, as.numeric(t(b) %*% solve(V) %*% b), tolerance = 1e-8)
})

test_that("with nothing imputed the degrees of freedom equal the design df exactly", {
  fit <- svyrcs(bmi ~ rcs(age, 4) + sex, design = mi_design_degenerate())
  expect_identical(unique(fit$curve$df), fit$degf)
  expect_identical(unique(fit$curve$fmi), 0)
  expect_identical(fit$tests$overall$df2, fit$degf)
  expect_identical(fit$tests$overall$fmi, 0)
})

test_that("an ordinary fit gains no df or fmi columns", {
  fit <- gaussian_fit()
  expect_false("df" %in% names(fit$curve))
  expect_false("fmi" %in% names(fit$curve))
  expect_null(fit$imputations)
  expect_null(fit$tests$overall$fmi)
})
