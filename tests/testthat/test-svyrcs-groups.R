test_that("svyrcs accepts an interaction and records the groups", {
  fit <- svyrcs(tchol ~ rcs(bmi, 4) * sex + age, design = nhanes_design())
  expect_s3_class(fit, "svyrcs")
  expect_equal(fit$groups$var, "sex")
  expect_equal(fit$groups$levels, c("Male", "Female"))
  expect_equal(fit$groups$ref_level, "Male")
  expect_equal(nrow(fit$curve), 200L * 2L)
  expect_true(all(is.finite(fit$curve$conf.low)))
})

test_that("grouped fits work for every model family", {
  design <- nhanes_design()
  skip_if_not_installed("survival")
  cases <- list(
    list(f = survival::Surv(time, event) ~ rcs(bmi, 4) * sex + age, family = NULL, measure = "HR"),
    list(f = tchol ~ rcs(bmi, 4) * sex + age, family = NULL, measure = "Difference"),
    list(f = high_chol ~ rcs(bmi, 4) * sex + age, family = quasibinomial(), measure = "OR"),
    list(f = statin_use ~ rcs(bmi, 4) * sex + age, family = quasipoisson(),
         measure = "Mean ratio")
  )
  for (cs in cases) {
    fit <- svyrcs(cs$f, design = design, family = cs$family)
    expect_equal(fit$measure, cs$measure, info = cs$measure)
    expect_equal(nlevels(fit$curve$group), 2L, info = cs$measure)
    expect_true(all(is.finite(fit$curve$estimate)), info = cs$measure)
    expect_true(is.finite(fit$tests$interaction$p_F), info = cs$measure)
  }
})

test_that("a four-level modifier gives four curves and nine interaction df", {
  fit <- svyrcs(tchol ~ rcs(bmi, 4) * race + age, design = nhanes_design(), n = 50)
  expect_equal(nlevels(fit$curve$group), 4L)
  expect_equal(nrow(fit$curve), 200L)
  expect_equal(fit$tests$interaction$df1, 9L)
  expect_length(fit$tests$by_group, 4L)
})

test_that("the reference is shared and every group passes through the null there", {
  fit <- svyrcs(tchol ~ rcs(bmi, 4) * race + age, design = nhanes_design())
  p <- predict(fit, x = fit$ref$value)
  expect_equal(p$estimate, rep(0, 4L))
  expect_equal(p$se, rep(0, 4L))
})

test_that("a subset design keeps its own degrees of freedom", {
  design <- nhanes_design()
  sub <- subset(design, cycle == "2005-2006")
  fit <- svyrcs(tchol ~ rcs(bmi, 4) * sex + age, design = sub)
  expect_equal(fit$degf, survey::degf(sub))
  expect_equal(fit$tests$interaction$df2, survey::degf(sub))
})

test_that("rows with a missing modifier are dropped and counted", {
  design <- nhanes_design()
  design$variables$grp <- nhanes_bmi$sex
  design$variables$grp[1:100] <- NA
  fit <- svyrcs(bmi ~ rcs(age, 4) * grp, design = design)
  expect_equal(fit$n_dropped, 100L)
  expect_equal(fit$n + fit$n_dropped, nrow(nhanes_bmi))
  expect_equal(fit$groups$var, "grp")
})

test_that("a rank-deficient interaction is diagnosed rather than reported as a missing column", {
  ## A modifier defined from the exposure makes the spline collinear within a group, so the fit
  ## drops interaction coefficients as aliased.
  design <- nhanes_design()
  design$variables$old <- factor(ifelse(nhanes_bmi$age >= 65, "old", "young"))
  expect_error(svyrcs(bmi ~ rcs(age, 4) * old, design = design), class = "svyrcs_error")
  expect_error(svyrcs(bmi ~ rcs(age, 4) * old, design = design), "rank deficient")
})

test_that("misuse of the interaction form is rejected", {
  design <- nhanes_design()
  expect_error(svyrcs(tchol ~ rcs(bmi, 4) * age, design = design), class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcs(bmi, 4) * sex * race, design = design), class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcs(bmi, 4):sex + sex, design = design), class = "svyrcs_error")
})

test_that("weighted knots are still shared across groups", {
  fit <- svyrcs(tchol ~ rcs(bmi, 4) * sex + age, design = nhanes_design())
  plain <- svyrcs(tchol ~ rcs(bmi, 4) + age, design = nhanes_design())
  expect_equal(fit$knots, plain$knots)
  expect_equal(fit$ref$value, plain$ref$value)
})
