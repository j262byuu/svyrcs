## What happens at the boundary with rms, recorded as fact rather than as a wish.
##
## Until 0.7.0 this package exported `rcs()`, which masked `rms::rcs()`. The rename removes the
## collision: `rcs()` in an rms formula is now unambiguously rms's. What remains is the converse --
## a user pasting `rcspline()` out of this package's vignette into an `ols()`, `lrm()` or `cph()`
## call. That is copy-paste, not a typo, so it will happen.
##
## The package does not currently prevent it. These tests pin what it does instead, so that if the
## behaviour ever changes -- either because rms changes, or because a future version adds a guard --
## somebody has to look at it deliberately.

skip_if_no_rms <- function() {
  skip_if_not_installed("rms")
  skip_if_not_installed("survival")
}

test_that("rcs() in an rms formula is rms's, now that we no longer export the name", {
  skip_if_no_rms()
  ## The whole point of the rename. Before it, whichever package was attached last won.
  expect_false("rcs" %in% getNamespaceExports("svyrcs"))
  expect_true("rcs" %in% getNamespaceExports("rms"))
})

test_that("a formula still written with rcs() is refused, and says where it would resolve", {
  des <- nhanes_design()
  expect_error(svyrcs(tchol ~ rcs(bmi, 4), design = des), class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcs(bmi, 4), design = des), "no longer exports")
  expect_error(svyrcs(tchol ~ rcs(bmi, 4), design = des), "rms::rcs")
  ## nested inside another term too, since that is where a rewrite would be easy to miss
  expect_error(svyrcs(tchol ~ age + rcs(bmi, 4), design = des), class = "svyrcs_error")
})

test_that("rcspline() inside an rms model fits, and loses the Nonlinear test", {
  ## RECORDING CURRENT BEHAVIOUR, NOT ENDORSING IT. rms identifies its spline terms through the
  ## Design attributes (`assume`, `nonlinear`, `parms`) that rms::rcs() attaches. This basis carries
  ## `knots`/`nk`/`var` instead, so rms treats it as a plain numeric matrix: the fit is correct, and
  ## the non-linearity test silently disappears.
  ##
  ## Whether to make this fail loudly is open -- see the v0.7 note in the task record. If it is ever
  ## addressed, this test is where the change surfaces.
  skip_if_no_rms()
  set.seed(2)
  n <- 500
  d <- data.frame(x = runif(n, 15, 45))
  d$y <- 0.02 * d$x - 0.0006 * (d$x - 30)^2 + rnorm(n, 0, 0.3)
  dd <- rms::datadist(d)
  old <- options(datadist = dd)
  on.exit(options(old), add = TRUE)

  m <- rms::ols(y ~ rcspline(x, 4), data = d)
  expect_s3_class(m, "ols")
  expect_length(stats::coef(m), 4L)

  ## the coefficient names lose the variable prefix, because rms builds it from the Design attributes
  expect_equal(unname(names(stats::coef(m))[-1L]), c("", "'", "''"))

  ## and anova() has no Nonlinear row
  a <- utils::capture.output(print(stats::anova(m)))
  expect_false(any(grepl("Nonlinear", a)))

  ## the fit itself is right: identical to the same model built with rms::rcs()
  m2 <- rms::ols(y ~ rms::rcs(x, 4), data = d)
  expect_equal(unname(stats::fitted(m)), unname(stats::fitted(m2)), tolerance = 1e-10)
  ## which is what makes the missing test the whole of the problem
  a2 <- utils::capture.output(print(stats::anova(m2)))
  expect_true(any(grepl("Nonlinear", a2)))
})

test_that("two rcspline() terms in an rms model are refused by rms's duplicate-name guard", {
  ## The basis columns are named "", "'", "''" with no variable prefix, so two spline terms collide.
  ## rms::Design() rejects the duplicate names, which is what stops the coefficients from being
  ## misattributed. Worth pinning: the safety here is upstream, not ours, and if rms ever
  ## disambiguates duplicates instead of erroring, silent misattribution becomes possible.
  skip_if_no_rms()
  set.seed(3)
  n <- 500
  d <- data.frame(x = runif(n, 15, 45), z = runif(n, 0, 10))
  d$y <- 0.02 * d$x + 0.15 * d$z + rnorm(n, 0, 0.3)
  dd <- rms::datadist(d)
  old <- options(datadist = dd)
  on.exit(options(old), add = TRUE)

  expect_error(rms::ols(y ~ rcspline(x, 4) + rcspline(z, 4), data = d), "duplicated column name")
})
