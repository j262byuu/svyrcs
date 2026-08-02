grouped_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- svyrcs(tchol ~ rcs(bmi, 4) * sex + age, design = nhanes_design())
    }
    cached
  }
})

test_that("the interaction F matches survey::regTermTest", {
  fit <- grouped_fit()
  int_term <- paste0(fit$term_label, ":", fit$groups$var)
  rt <- survey::regTermTest(fit$model, reformulate(int_term), df = fit$degf)
  expect_equal(as.numeric(rt$Ftest), fit$tests$interaction$F, tolerance = 1e-8)
  expect_equal(as.numeric(rt$df), fit$tests$interaction$df1)
  expect_equal(as.numeric(rt$p), fit$tests$interaction$p_F, tolerance = 1e-8)
})

test_that("the interaction and shape tests have the right degrees of freedom", {
  fit <- grouped_fit()
  nk <- fit$nk
  nlev <- length(fit$groups$levels)
  expect_equal(fit$tests$interaction$df1, (nk - 1L) * (nlev - 1L))
  expect_equal(fit$tests$shape$df1, (nk - 2L) * (nlev - 1L))
  expect_equal(fit$tests$interaction$df2, fit$degf)

  four <- svyrcs(tchol ~ rcs(bmi, 4) * race + age, design = nhanes_design())
  expect_equal(four$tests$interaction$df1, 3L * 3L)
  expect_equal(four$tests$shape$df1, 2L * 3L)
})

test_that("the joint overall test spans the spline block and its interactions", {
  fit <- grouped_fit()
  nk <- fit$nk
  nlev <- length(fit$groups$levels)
  expect_equal(fit$tests$overall$df1, (nk - 1L) * nlev)
  expect_equal(fit$tests$nonlinear$df1, (nk - 2L) * nlev)
})

test_that("per-group tests match a manual Wald test on that group's coefficients", {
  fit <- grouped_fit()
  term <- find_rcs_term(fit$model)
  m <- find_modifier(fit$model, term)
  b <- coef(fit$model)
  V <- vcov(fit$model)

  for (g in fit$groups$levels) {
    L <- selection_matrix(names(b), m$spline_cols, m$columns[[g]])
    bg <- as.numeric(L %*% b)
    Vg <- L %*% V %*% t(L)
    chisq <- as.numeric(t(bg) %*% solve(Vg) %*% bg)
    expect_equal(fit$tests$by_group[[g]]$overall$chisq, chisq, tolerance = 1e-8, info = g)
  }
})

test_that("the reference group's tests use the main effects alone", {
  fit <- grouped_fit()
  term <- find_rcs_term(fit$model)
  m <- find_modifier(fit$model, term)
  b <- coef(fit$model)[m$spline_cols]
  V <- vcov(fit$model)[m$spline_cols, m$spline_cols, drop = FALSE]
  manual <- rcs_tests(b, V, fit$degf)
  expect_equal(fit$tests$by_group[[fit$groups$ref_level]]$overall$chisq,
               manual$overall$chisq, tolerance = 1e-10)
})

test_that("an ungrouped fit keeps the 0.1.0 tests shape", {
  fit <- gaussian_fit()
  expect_named(fit$tests, c("overall", "nonlinear"))
  expect_null(fit$tests$interaction)
  expect_null(fit$tests$by_group)
})

test_that("both formula orders give the same interaction test", {
  a <- svyrcs(tchol ~ rcs(bmi, 4) * sex + age, design = nhanes_design())
  b <- svyrcs(tchol ~ sex * rcs(bmi, 4) + age, design = nhanes_design())
  expect_equal(a$tests$interaction$F, b$tests$interaction$F, tolerance = 1e-10)
  expect_equal(a$tests$shape$F, b$tests$shape$F, tolerance = 1e-10)
  expect_equal(as.data.frame(a$curve), as.data.frame(b$curve))
})
