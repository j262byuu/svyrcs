test_that("all four model families work on an imputed design", {
  skip_if_not_installed("survival")
  designs <- mi_design()
  cases <- list(
    list(f = survival::Surv(time, event) ~ rcs(bmi, 4) + age + tchol, family = NULL,
         measure = "HR"),
    list(f = bmi ~ rcs(age, 4) + sex + tchol, family = NULL, measure = "Difference"),
    list(f = high_chol ~ rcs(bmi, 4) + age + hdl, family = quasibinomial(), measure = "OR"),
    list(f = statin_use ~ rcs(bmi, 4) + age + tchol, family = quasipoisson(), measure = "RR")
  )
  for (cs in cases) {
    fit <- svyrcs(cs$f, design = designs, family = cs$family)
    expect_equal(fit$measure, cs$measure, info = cs$measure)
    expect_equal(fit$imputations$m, 5L, info = cs$measure)
    expect_true(all(is.finite(fit$curve$estimate)), info = cs$measure)
    expect_true(all(is.finite(fit$curve$conf.high)), info = cs$measure)
    expect_true(all(fit$curve$df <= fit$degf + 1e-9), info = cs$measure)
  }
})

test_that("the interaction form works under imputation", {
  fit <- svyrcs(bmi ~ rcs(age, 4) * sex + tchol, design = mi_design())
  expect_equal(fit$groups$var, "sex")
  expect_equal(nrow(fit$curve), 400L)
  expect_true(all(c("df", "fmi") %in% names(fit$curve)))
  expect_lte(fit$tests$interaction$df2, fit$degf + 1e-9)
  expect_lte(fit$tests$shape$df2, fit$degf + 1e-9)
  expect_length(fit$tests$by_group, 2L)
  for (g in fit$groups$levels) {
    expect_lte(fit$tests$by_group[[g]]$overall$df2, fit$degf + 1e-9)
  }
})

test_that("knots are the average of the per-imputation weighted quantiles", {
  designs <- mi_design()
  fit <- svyrcs(bmi ~ rcs(tchol, 4) + age, design = designs)

  each <- vapply(designs$designs, function(d) {
    as.numeric(coef(survey::svyquantile(~tchol, d, c(0.05, 0.35, 0.65, 0.95),
                                        na.rm = TRUE, ci = FALSE)))
  }, numeric(4))
  expect_equal(fit$knots, signif(rowMeans(each), 7))

  ## and every component fit used those same knots
  for (f in fit$model$fits) {
    expect_equal(find_rcs_term(f)$knots, fit$knots)
  }
})

test_that("the reference is pooled across imputations", {
  designs <- mi_design()
  fit <- svyrcs(bmi ~ rcs(tchol, 4) + age, design = designs)
  each <- vapply(designs$designs, function(d) {
    as.numeric(coef(survey::svyquantile(~tchol, d, 0.5, na.rm = TRUE, ci = FALSE)))
  }, numeric(1))
  expect_equal(fit$ref$value, mean(each))
  expect_gt(diff(range(each)), 0)  # the imputations really do disagree
})

test_that("with nothing imputed the fit reproduces the complete-data fit exactly", {
  plain <- svyrcs(bmi ~ rcs(age, 4) + sex, design = nhanes_design())
  imputed <- svyrcs(bmi ~ rcs(age, 4) + sex, design = mi_design_degenerate())

  expect_equal(imputed$knots, plain$knots)
  expect_equal(imputed$ref$value, plain$ref$value)
  expect_equal(imputed$n, plain$n)
  expect_equal(imputed$degf, plain$degf)
  expect_equal(coef(imputed$model), coef(plain$model))
  expect_equal(unname(vcov(imputed$model)), unname(vcov(plain$model)))

  expect_equal(imputed$curve$estimate, plain$curve$estimate)
  expect_equal(imputed$curve$conf.low, plain$curve$conf.low)
  expect_equal(imputed$curve$conf.high, plain$curve$conf.high)
  expect_equal(imputed$tests$overall$F, plain$tests$overall$F)
  expect_equal(imputed$tests$overall$p_F, plain$tests$overall$p_F)
  expect_equal(imputed$tests$nonlinear$p_F, plain$tests$nonlinear$p_F)
})

test_that("predict works on an imputed fit and carries the df through", {
  fit <- svyrcs(bmi ~ rcs(age, 4) + sex + tchol, design = mi_design())
  p <- predict(fit, x = c(30, 50, 70))
  expect_equal(nrow(p), 3L)
  expect_true(all(c("df", "fmi") %in% names(p)))
  expect_true(all(p$df <= fit$degf + 1e-9))
  expect_equal(predict(fit, x = fit$ref$value)$estimate, 0)
})

test_that("a subset of the imputed design keeps its own degrees of freedom", {
  sub <- subset(mi_design(), cycle == "2005-2006")
  fit <- svyrcs(bmi ~ rcs(age, 4) + sex + tchol, design = sub)
  expect_equal(fit$degf, survey::degf(sub$designs[[1]]))
  expect_lt(fit$degf, survey::degf(mi_design()$designs[[1]]))
})

test_that("malformed imputed designs are rejected", {
  one <- survey::svydesign(id = ~psu, strata = ~strata, weights = ~weight, nest = TRUE,
                           data = mitools::imputationList(list(nhanes_bmi)))
  expect_error(svyrcs(bmi ~ rcs(age, 4), design = one), class = "svyrcs_error")
  expect_error(svyrcs(bmi ~ rcs(age, 4), design = mitools::imputationList(list(nhanes_bmi))),
               class = "svyrcs_error")
  expect_error(svyrcs(bmi ~ rcs(nosuchvar, 4), design = mi_design()), class = "svyrcs_error")
})

test_that("the pooled model is inspectable but refuses to predict", {
  fit <- svyrcs(bmi ~ rcs(age, 4) + sex + tchol, design = mi_design())
  expect_s3_class(fit$model, "svyrcs_mifit")
  expect_length(fit$model$fits, 5L)
  expect_s3_class(fit$model$fits[[1]], "svyglm")
  expect_error(predict(fit$model, newdata = nhanes_bmi[1:2, ]), class = "svyrcs_error")
  expect_output(print(fit$model), "Pooled fit across 5 imputations")
})

test_that("print reports the imputation count and the fraction of missing information", {
  fit <- svyrcs(bmi ~ rcs(age, 4) + sex + tchol, design = mi_design())
  out <- capture.output(print(fit))
  expect_true(any(grepl("Imputed    m = 5", out)))
  expect_true(any(grepl("fraction of missing information", out)))
  expect_true(any(grepl("Barnard-Rubin", out)))
  ## a non-integer denominator df is shown to one decimal rather than rounded away
  expect_true(any(grepl("and [0-9]+\\.[0-9] df", out)))
})

test_that("an ordinary fit says nothing about imputation", {
  out <- capture.output(print(gaussian_fit()))
  expect_false(any(grepl("Imputed", out)))
  expect_false(any(grepl("Barnard-Rubin", out)))
})

test_that("summary and plot work on an imputed fit", {
  fit <- svyrcs(bmi ~ rcs(age, 4) + sex + tchol, design = mi_design())
  expect_output(print(summary(fit)), "Curve shape")
  p <- autoplot(fit)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 200L)
  expect_silent(invisible(ggplot2::ggplot_build(p)))
})
