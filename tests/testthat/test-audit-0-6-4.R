## Regressions for the 0.6.3 audit findings. Each is built from the reproduction that exposed it;
## the audit record is in .trellis/tasks/archive/2026-08/08-05-codex-audit/research/findings.md.

test_that("a Cox model fits on a replicate-weight design", {
  ## The package advertises replicate designs and Cox models; their intersection had never worked.
  ## survey's svycoxph.svyrep.design records match.call() and evaluates it as with(data, eval(g)),
  ## so a formula passed by variable name is out of scope by the time it is used.
  skip_if_not_installed("survival")
  rd <- rep_design()
  fit <- suppressWarnings(
    svyrcs(survival::Surv(time, event) ~ rcs(bmi, 4) + age, design = rd, n = 5))
  expect_s3_class(fit, "svyrcs")
  expect_true(all(is.finite(coef(fit$model))))

  ## and it agrees with the linearised fit on the same knots, which is the point of the fix
  lin <- suppressWarnings(
    svyrcs(survival::Surv(time, event) ~ rcs(bmi, 4) + age, design = nhanes_design(), n = 5))
  expect_equal(unname(coef(fit$model)), unname(coef(lin$model)), tolerance = 1e-10)

  ## the GLM path on the same design was already fine and must stay so
  expect_s3_class(suppressWarnings(svyrcs(tchol ~ rcs(bmi, 4) + age, design = rd, n = 5)), "svyrcs")
})

test_that("a knot count that cannot be placed is refused, not silently shortened", {
  ## small_sample_outer_knots() ended in sort(unique(kn)), so when the fifth-order-statistic
  ## replacements collided with the interior knots the vector shrank and the caller fitted a
  ## different spline than it asked for.
  expect_error(rcs_knots(1:9, 5), class = "svyrcs_error")
  expect_error(rcs_knots(1:9, 5), "coincide with the interior knots")
  expect_error(rcs_knots(1:3, 3), class = "svyrcs_error")
  expect_error(rcs_knots(1:4, 4), class = "svyrcs_error")
  expect_error(rcs_knots(1:4, 4), "at least 5")

  ## n = 8 crosses the two indices (5 and 4) but still yields four distinct knots, and Hmisc agrees,
  ## so the guard is set by where the rule is undefined rather than by where it looks tidy
  expect_length(rcs_knots(1:8, 4), 4L)
  expect_equal(rcs_knots(1:8, 4), c(3.45, 4, 5, 5.55))
  expect_length(rcs_knots(1:20, 4), 4L)
})

test_that("a non-syntactic exposure name works", {
  ## stats::reformulate() runs str2lang() on the bare label, so `waist cm` was a parse error -- on
  ## the default path, since both weighted helpers went through it.
  d <- nhanes_bmi[stats::complete.cases(
    nhanes_bmi[c("bmi", "tchol", "psu", "strata", "weight")]), ]
  d[["waist cm"]] <- d$bmi * 2.5
  des <- survey::svydesign(id = ~psu, strata = ~strata, weights = ~weight, nest = TRUE, data = d)

  fit <- svyrcs(tchol ~ rcs(`waist cm`, 4), design = des, n = 5)
  expect_length(fit$knots, 4L)
  expect_true(all(is.finite(fit$knots)))

  ## the weighted helpers specifically, since that is where it failed
  expect_true(is.finite(svyrcs:::weighted_mean("waist cm", des)))
  expect_length(svyrcs:::weighted_quantile("waist cm", des, c(0.25, 0.75)), 2L)
})

test_that("a numeric NA level is refused with a typed error", {
  ## The 0.6.2 test pinned `level = NA`, which is logical and fails is.numeric(). NA_real_ and NaN
  ## took a different branch and reached the base "missing value where TRUE/FALSE needed".
  des <- nhanes_design()
  fit <- svyrcs(tchol ~ rcs(bmi, 4) + age, design = des, n = 5)
  for (bad in list(NA, NA_real_, NaN, Inf, -Inf)) {
    expect_error(predict(fit, x = 30, level = bad), class = "svyrcs_error")
  }
  expect_error(predict(fit, x = 30, level = NA_real_), "finite")
})

test_that("`at` and `range` are not both accepted", {
  ## `at` was documented as making `range` ignored, and it was -- for the grid but not for the
  ## reference search, so at = c(20, 25, 30) with range = c(100, 200) and ref = "min" returned a
  ## grid of 20-30 anchored at 200.
  des <- nhanes_design()
  fit <- svyrcs(tchol ~ rcs(bmi, 4) + age, design = des, n = 5)
  expect_error(svyrcs_curve(fit$model, var = "bmi", at = c(20, 25, 30), range = c(100, 200)),
               class = "svyrcs_error")
  expect_error(svyrcs_curve(fit$model, var = "bmi", at = c(20, 25, 30), range = c(100, 200)),
               "not both")

  ## either alone still works, and `at` alone anchors inside its own values
  cv <- svyrcs_curve(fit$model, var = "bmi", at = c(20, 25, 30), ref = "min")
  expect_equal(range(cv$x), c(20, 30))
  expect_true(attr(cv, "ref") >= 20 && attr(cv, "ref") <= 30)
  expect_equal(nrow(svyrcs_curve(fit$model, var = "bmi", range = c(20, 40), n = 6)), 6L)
})

test_that("an invalid imputation denominator is refused rather than turned into NaN", {
  ## Reiter's expansion returns a negative denominator once the between-imputation variance
  ## dominates -- v crosses zero around FMI 0.893 -- and pf() then gave NaN with only base R's
  ## "NaNs produced" to show for it. No construction reached this from a proper imputation, so this
  ## guards an arithmetic possibility.
  k <- 3
  Q <- rep(0.4, k); Ubar <- diag(k) * 0.1

  good <- svyrcs:::d1_test(Q, Ubar, diag(k) * 0.1, 20, dfcom = 31)
  expect_true(is.finite(good$v) && good$v > 0)

  bad <- svyrcs:::d1_test(Q, Ubar, diag(k) * 0.1 * 1e6, 20, dfcom = 31)
  expect_s3_class(bad, "svyrcs_invalid_df")
  expect_lt(bad$v, 0)

  ## and it surfaces as an unavailable test with its own reason, not as a NaN p-value
  mi <- list(Ubar = Ubar, B = diag(k) * 0.1 * 1e6, m = 20, degf = 31)
  tt <- suppressWarnings(svyrcs:::wald_block(Q, Ubar + mi$B, 31, "overall association", mi))
  expect_equal(tt$reason, "invalid-df")
  expect_true(is.na(tt$p_F))
  expect_warning(svyrcs:::wald_block(Q, Ubar + mi$B, 31, "overall association", mi),
                 "not a valid reference distribution")
  expect_match(svyrcs:::fmt_test(tt, "Overall association"), "invalid denominator")
})

test_that("the small-sample knot rule says when it does not apply", {
  ## Harrell's rule is defined on order statistics of an unweighted sample. The weighted path -- the
  ## default -- cannot apply it, and 0.6.2 claimed the default matched rms below n = 100 on the
  ## strength of tests that all went through the unweighted path.
  set.seed(1)
  n <- 60
  d <- data.frame(x = sort(rgamma(n, 2, 0.4)), psu = rep(1:20, each = 3),
                  strata = rep(1:10, each = 6), w = 1)
  d$y <- 0.3 * d$x + rnorm(n)
  des <- survey::svydesign(id = ~psu, strata = ~strata, weights = ~w, nest = TRUE, data = d)

  expect_warning(fw <- svyrcs(y ~ rcs(x, 4), design = des, n = 5), "60 analytic observations")
  expect_warning(svyrcs(y ~ rcs(x, 4), design = des, n = 5), "weighted_knots = FALSE")

  ## the unweighted path applies the rule and does not warn about it
  fu <- svyrcs(y ~ rcs(x, 4), design = des, n = 5, weighted_knots = FALSE)
  expect_equal(fu$knots[c(1, 4)], c(sort(d$x)[5], sort(d$x)[n - 4]), tolerance = 1e-6)
  expect_false(isTRUE(all.equal(fw$knots, fu$knots)))

  ## and a large sample warns about neither
  expect_silent(svyrcs(tchol ~ rcs(bmi, 4), design = nhanes_design(), n = 5))
})
