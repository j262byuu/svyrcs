test_that("the Cox path gives hazard ratios and counts observations, not events", {
  skip_if_not_installed("survival")
  fit <- cox_fit()
  expect_s3_class(fit, "svyrcs")
  expect_s3_class(fit$model, "svycoxph")
  expect_equal(fit$measure, "HR")
  expect_equal(fit$null, 1)
  ## nobs() on a coxph object returns the event count, which would be badly wrong here
  expect_equal(fit$n, nrow(nhanes_bmi))
  expect_equal(fit$nevents, sum(nhanes_bmi$event))
  expect_lt(fit$nevents, fit$n)
  expect_equal(fit$degf, survey::degf(nhanes_design()))
})

test_that("each glm family maps to the right effect measure", {
  design <- nhanes_design()
  cases <- list(
    list(f = tchol ~ rcspline(bmi, 4) + age, family = NULL, measure = "Difference", null = 0),
    list(f = high_chol ~ rcspline(bmi, 4) + age, family = quasibinomial(), measure = "OR", null = 1),
    list(f = statin_use ~ rcspline(bmi, 4) + age, family = quasipoisson(), measure = "Mean ratio",
         null = 1)
  )
  for (cs in cases) {
    fit <- svyrcs(cs$f, design = design, family = cs$family)
    expect_equal(fit$measure, cs$measure, info = cs$measure)
    expect_equal(fit$null, cs$null, info = cs$measure)
    expect_true(all(is.finite(fit$curve$estimate)), info = cs$measure)
    expect_true(all(is.finite(fit$curve$conf.low)), info = cs$measure)
    expect_true(all(is.finite(fit$curve$conf.high)), info = cs$measure)
    expect_true(is.finite(fit$tests$overall$p_F), info = cs$measure)
  }
})

test_that("weighted knots differ from unweighted ones and both are reproducible", {
  design <- nhanes_design()
  fw <- svyrcs(tchol ~ rcspline(bmi, 4) + age, design = design, weighted_knots = TRUE)
  fu <- svyrcs(tchol ~ rcspline(bmi, 4) + age, design = design, weighted_knots = FALSE)

  ## Both are computed on the analytic sample -- the rows the model is fitted on -- not on the
  ## design as supplied; `tchol` is missing for 649 of them.
  keep <- stats::complete.cases(nhanes_bmi[c("tchol", "bmi", "age")])
  analytic <- design[keep, ]
  expect_equal(fu$knots,
               signif(unname(quantile(nhanes_bmi$bmi[keep], c(0.05, 0.35, 0.65, 0.95))), 7))
  expect_equal(fw$knots,
               signif(unname(coef(survey::svyquantile(~bmi, analytic, c(0.05, 0.35, 0.65, 0.95),
                                                      na.rm = TRUE, ci = FALSE))), 7))
  expect_false(isTRUE(all.equal(fw$knots, fu$knots)))
  expect_true(fw$weighted_knots)
  expect_false(fu$weighted_knots)
})

test_that("explicit knots are used verbatim regardless of weighted_knots", {
  kn <- c(20, 25, 30, 40)
  for (w in c(TRUE, FALSE)) {
    fit <- svyrcs(tchol ~ rcspline(bmi, kn) + age, design = nhanes_design(), weighted_knots = w)
    expect_equal(fit$knots, kn)
    expect_equal(fit$nk, 4L)
  }
})

test_that("the fitted model carries explicit knots and is self-contained", {
  fit <- gaussian_fit()
  ## no symbol lookup needed: the knots are inlined in the terms, so they evaluate with nothing
  ## but base R in scope
  pv <- as.list(attr(terms(fit$model), "predvars"))
  is_rcs <- vapply(pv, function(e) is.call(e) &&
                     identical(as.character(e[[1L]])[1L], "rcspline"), logical(1L))
  expect_true(any(is_rcs))
  expect_equal(eval(pv[is_rcs][[1L]]$knots, envir = baseenv()), fit$knots)

  reloaded <- suppressWarnings(unserialize(serialize(fit$model, NULL)))
  expect_equal(coef(reloaded), coef(fit$model))
  nd <- nhanes_bmi[1:3, ]
  nd$bmi <- c(18, 28, 45)
  expect_equal(as.numeric(predict(reloaded, newdata = nd, type = "link")),
               as.numeric(predict(fit$model, newdata = nd, type = "link")))
})

test_that("the number of knots changes the model as expected", {
  design <- nhanes_design()
  for (nk in 3:7) {
    fit <- svyrcs(tchol ~ rcspline(bmi, nk) + age, design = design)
    expect_equal(fit$nk, nk)
    expect_length(fit$knots, nk)
    expect_equal(fit$tests$overall$df1, nk - 1L)
    expect_length(spline_coef_index(fit$model, find_rcs_term(fit$model)), nk - 1L)
  }
})

test_that("a subset design keeps its own degrees of freedom and sample size", {
  design <- nhanes_design()
  sub <- subset(design, cycle == "2005-2006")
  fit <- svyrcs(tchol ~ rcspline(bmi, 4) + age, design = sub)
  expect_equal(fit$degf, survey::degf(sub))
  expect_lt(fit$degf, survey::degf(design))
  expect_lt(fit$n, nrow(nhanes_bmi))
})

test_that("rows dropped for missing data are counted and reported", {
  fit <- gaussian_fit()
  expect_equal(fit$n_dropped, sum(is.na(nhanes_bmi$tchol)))
  expect_equal(fit$n + fit$n_dropped, nrow(nhanes_bmi))
  expect_output(print(fit), "Dropped")

  cox <- cox_fit()
  expect_equal(cox$n_dropped, 0L)
  expect_false(any(grepl("Dropped", capture.output(print(cox)))))
})

test_that("misuse is rejected with informative errors", {
  design <- nhanes_design()
  expect_error(svyrcs(tchol ~ age + sex, design = design), class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4) + rcspline(age, 4), design = design), class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = nhanes_bmi), class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(nosuchvar, 4), design = design), class = "svyrcs_error")
  expect_error(svyrcs(~ rcspline(bmi, 4), design = design), class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(bmi, 2), design = design), class = "svyrcs_error")
})

test_that("family is ignored with a warning for a Cox model", {
  skip_if_not_installed("survival")
  expect_warning(
    svyrcs(survival::Surv(time, event) ~ rcspline(bmi, 4) + age, design = nhanes_design(),
           family = quasibinomial()),
    "ignored"
  )
})

test_that("a bare Surv() response works when survival is attached", {
  skip_if_not_installed("survival")
  library(survival)
  fit <- svyrcs(Surv(time, event) ~ rcspline(bmi, 4) + age, design = nhanes_design())
  expect_s3_class(fit$model, "svycoxph")
  expect_equal(fit$measure, "HR")
})
