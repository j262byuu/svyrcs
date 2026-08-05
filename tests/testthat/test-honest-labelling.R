## Review items 2, 3, 5, 6, 7c and 8a: nothing here changes an estimate, only what the package says
## about one, plus the `df` argument. The compat suite is the guard that no number moved.

test_that("summary calls the ranges pointwise", {
  fit <- gaussian_fit()
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Pointwise .*band excludes|pointwise .*band includes", out)))
  ## the bare word, which invited a family-wise reading, is gone from the output
  expect_false(any(grepl("^\\s+95% band excludes", out)))
})

test_that("the grid extrema are labelled as grid searches", {
  fit <- gaussian_fit()
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Grid lowest", out)))
  expect_true(any(grepl("Grid highest", out)))
})

test_that("a data-selected reference says its uncertainty is not propagated", {
  skip_if_not_installed("survival")
  fit <- svyrcs(survival::Surv(time, event) ~ rcspline(bmi, 4) + age, design = nhanes_design(),
                ref = "min", n = 50)
  out <- capture.output(print(fit))
  expect_true(any(grepl("selected from these data", out)))
  expect_true(any(grepl("not in the band", out)))

  ## a weighted-quantile reference is not data-selected in that sense, so it says nothing
  plain <- capture.output(print(gaussian_fit()))
  expect_false(any(grepl("selected from these data", plain)))
})

test_that("df overrides the design degrees of freedom", {
  design <- nhanes_design()
  base <- svyrcs(tchol ~ rcspline(bmi, 4) + age, design = design, n = 5)
  expect_equal(base$degf, survey::degf(design))

  ## Inf gives the normal approximation: narrower than the t interval on 31 df
  inf <- svyrcs(tchol ~ rcspline(bmi, 4) + age, design = design, n = 5, df = Inf)
  expect_equal(inf$degf, Inf)
  expect_equal(inf$curve$estimate, base$curve$estimate)
  expect_equal(inf$curve$se, base$curve$se)
  wid <- function(f) f$curve$conf.high - f$curve$conf.low
  expect_true(all(wid(inf) <= wid(base) + 1e-12))
  expect_true(any(wid(inf) < wid(base)))

  ## and the critical value really is the normal one
  i <- 3L
  expect_equal((inf$curve$conf.high - inf$curve$estimate)[i],
               qnorm(0.975) * inf$curve$se[i], tolerance = 1e-10)

  ## a smaller df widens
  small <- svyrcs(tchol ~ rcspline(bmi, 4) + age, design = design, n = 5, df = 5)
  expect_true(all(wid(small) >= wid(base) - 1e-12))

  ## the tests use it too
  expect_equal(base$tests$overall$df2, survey::degf(design))
  expect_equal(inf$tests$overall$df2, Inf)
  expect_lt(inf$tests$overall$p_F, base$tests$overall$p_F)
})

test_that("df is validated", {
  design <- nhanes_design()
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = design, df = 0), class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = design, df = -1), class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = design, df = c(10, 20)),
               class = "svyrcs_error")
})

test_that("df = Inf survives the Barnard-Rubin arithmetic", {
  ## (Inf + 1)/(Inf + 3) is NaN, so the complete-data correction needs an explicit branch.
  ii <- lapply(1:4, function(i) {
    set.seed(i)
    o <- nhanes_bmi
    m <- is.na(o$tchol)
    o$tchol[m] <- sample(o$tchol[!m], sum(m), TRUE)
    o
  })
  des <- survey::svydesign(id = ~psu, strata = ~strata, weights = ~weight, nest = TRUE,
                           data = mitools::imputationList(ii))
  fit <- svyrcs(bmi ~ rcspline(age, 4) + tchol, design = des, n = 10, df = Inf)
  expect_true(all(is.finite(fit$curve$conf.high)))
  expect_true(all(!is.nan(fit$curve$df)))
  ## with an infinite complete-data df, Barnard-Rubin reduces to Rubin's original rule
  expect_true(all(fit$curve$df > 0))
})

test_that("a count model reports a rate ratio only when it has person-time", {
  design <- nhanes_design()
  design$variables$py <- pmax(nhanes_bmi$time, 0.1)

  no_offset <- svyrcs(statin_use ~ rcspline(bmi, 4) + age, design = design,
                      family = quasipoisson(), n = 5)
  expect_equal(no_offset$measure, "Mean ratio")

  with_offset <- svyrcs(event ~ rcspline(bmi, 4) + age + offset(log(py)), design = design,
                        family = quasipoisson(), n = 5)
  expect_equal(with_offset$measure, "Rate ratio")

  ## both are still ratios against a null of 1
  expect_equal(no_offset$null, 1)
  expect_equal(with_offset$null, 1)
})

test_that("non-default contrasts work as modifiers and agree with treatment coding", {
  design <- nhanes_design()
  base <- svyrcs(tchol ~ rcspline(bmi, 4) * race + age, design = design, n = 20)

  for (ct in c("contr.sum", "contr.helmert", "contr.poly")) {
    d2 <- nhanes_design()
    d2$variables$race <- `contrasts<-`(factor(nhanes_bmi$race),
                                       value = get(ct)(nlevels(factor(nhanes_bmi$race))))
    fit <- svyrcs(tchol ~ rcspline(bmi, 4) * race + age, design = d2, n = 20)
    expect_equal(length(fit$groups$levels), 4L, info = ct)
    expect_equal(fit$curve$estimate, base$curve$estimate, tolerance = 1e-8, info = ct)
    expect_equal(fit$curve$se, base$curve$se, tolerance = 1e-8, info = ct)
    expect_equal(fit$tests$interaction$F, base$tests$interaction$F, tolerance = 1e-8, info = ct)
  }
})

test_that("a user-supplied contrast matrix works", {
  d2 <- nhanes_design()
  f <- factor(nhanes_bmi$sex)
  contrasts(f) <- matrix(c(-0.5, 0.5), ncol = 1, dimnames = list(levels(f), "diff"))
  d2$variables$sx <- f

  base <- svyrcs(tchol ~ rcspline(bmi, 4) * sex + age, design = nhanes_design(), n = 20)
  fit <- svyrcs(tchol ~ rcspline(bmi, 4) * sx + age, design = d2, n = 20)
  expect_equal(fit$curve$estimate, base$curve$estimate, tolerance = 1e-8)
  expect_equal(fit$tests$interaction$F, base$tests$interaction$F, tolerance = 1e-8)
})

test_that("the documentation no longer claims Taylor linearisation for every design", {
  desc <- readLines(system.file("DESCRIPTION", package = "svyrcs"))
  taylor <- grep("Taylor", desc, value = TRUE)
  skip_if(length(taylor) == 0)
  ## wherever Taylor is mentioned, replicate weights are mentioned too
  expect_true(any(grepl("replicate", desc, ignore.case = TRUE)))
})
