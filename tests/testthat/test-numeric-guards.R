## Guards for the three numerical concerns that were left open after the 0.6.0 review and settled by
## adversarial construction in 0.6.3. Each test is built from the script that produced the verdict.

test_that("nearly coincident knots are named as such", {
  ## anyDuplicated() catches exact ties. What slips through is a cluster with floating-point jitter:
  ## a biomarker below the limit of detection substituted by a constant, where the substitution
  ## arithmetic leaves the values distinct around 1e-15.
  set.seed(21)
  n <- 4000
  lod <- 0.35
  below <- n * 0.7
  x <- c(rep(lod / sqrt(2), below) * (1 + rnorm(below, 0, 1e-15)),
         lod + rgamma(n - below, 2, 0.6))

  expect_warning(kn <- rcs_knots(x, 4), "almost coincident")
  ## distinct, so the exact-duplicate check cannot see them
  expect_false(anyDuplicated(kn) > 0)
  expect_lt(min(diff(kn)) / diff(range(kn)), 1e-6)

  ## the message names the spacing rather than blaming the design
  expect_warning(rcs_knots(x, 4), "knot range")
  expect_warning(rcs_knots(x, 4), "limit of detection")
})

test_that("explicit knots get the same spacing check as derived ones", {
  expect_warning(rcs_knots(1:100, c(10, 45, 45 + 1e-9, 90)), "almost coincident")
  expect_silent(rcs_knots(1:100, c(10, 45, 60, 90)))
})

test_that("the spacing warning stays quiet on ordinary exposures", {
  ## The floor is the false-alarm sweep: a warning that fires on real data is one users learn to
  ## ignore, and then it guards nothing.
  set.seed(5)
  n <- 2000
  shapes <- list(runif(n, 0, 100), rnorm(n, 50, 10), rgamma(n, 2, 0.5), rgamma(n, 1, 0.3),
                 rlnorm(n, 0, 1), rexp(n, 1), nhanes_bmi$bmi[!is.na(nhanes_bmi$bmi)])
  for (x in shapes) {
    for (nk in c(4, 5, 6)) {
      expect_silent(rcs_knots(x, nk))
    }
  }
})

test_that("a non-finite variance is reported as a variance problem, not a convergence problem", {
  ## Reachable through a replicate design whose rscales carry negative coefficients: survey takes
  ## their square root, so vcov() comes back with NaN while the coefficients stay finite. Before
  ## 0.6.3 the only sign was base R's "NaNs produced" from sqrt().
  skip_if_not_installed("survey")
  des <- nan_vcov_design()

  w <- capture_warnings(fit <- svyrcs(tchol ~ rcs(bmi, 5) + age, design = des, n = 20))

  expect_true(all(is.finite(fit$curve$estimate)))
  expect_true(all(is.na(fit$curve$se)))
  expect_match(w, "non-finite standard errors", all = FALSE)
  ## and it must not blame the model for failing to converge, which it did not
  expect_false(any(grepl("converge to finite coefficients", w)))
})

test_that("an unavailable test says which of its causes applies", {
  skip_if_not_installed("survey")
  des <- nan_vcov_design()
  suppressWarnings(fit <- svyrcs(tchol ~ rcs(bmi, 5) + age, design = des, n = 5))

  ## four spline columns are in the model, so "a linear term only" would be plainly false
  expect_equal(length(grep("rcs\\(", names(coef(fit$model)))), 4L)
  expect_false(all(is.finite(vcov(fit$model))))

  out <- capture.output(print(fit))
  expect_match(paste(out, collapse = "\n"), "not estimable \\(the covariance block is not finite\\)")
  expect_false(any(grepl("a linear term only", out)))

  expect_equal(fit$tests$overall$reason, "non-finite")
})

test_that("a genuine linear-term-only block still says so", {
  des <- nhanes_design()
  ## three knots leaves two spline columns; dropping the linear one leaves a single column, so the
  ## non-linearity block is still estimable. The linear reason is reachable from rcs_tests() when the
  ## block has no columns at all.
  tt <- svyrcs:::rcs_tests(c(x = 1), matrix(1, 1, 1), 30)
  expect_equal(tt$nonlinear$reason, "linear")
  expect_match(svyrcs:::fmt_test(tt$nonlinear, "Non-linearity"), "a linear term only")
})

test_that("every wald_block path carries a reason", {
  ok <- svyrcs:::wald_block(c(1, 2), diag(2), 30, "x")
  expect_true(is.na(ok$reason))

  nf <- suppressWarnings(svyrcs:::wald_block(c(1, 2), matrix(c(1, NaN, NaN, 1), 2), 30, "x"))
  expect_equal(nf$reason, "non-finite")
  expect_warning(svyrcs:::wald_block(c(1, 2), matrix(c(1, NaN, NaN, 1), 2), 30, "x"),
                 "non-finite values")

  sing <- suppressWarnings(svyrcs:::wald_block(c(1, 2), matrix(1, 2, 2), 30, "x"))
  expect_equal(sing$reason, "singular")
})

test_that("the rank guard fires before the Wald arithmetic degrades", {
  ## C2 was refuted because warn_if_rank_deficient() gives way at a far lower condition number than
  ## solve() does. That margin is a property of qr()'s default tolerance and nothing recorded the
  ## dependency, so a future change to it could open the gap this test exists to keep shut.
  make_cov <- function(k, cond, seed = 11) {
    set.seed(seed)
    Q <- qr.Q(qr(matrix(rnorm(k * k), k, k)))
    V <- Q %*% diag(exp(seq(0, -log(cond), length.out = k))) %*% t(Q)
    (V + t(V)) / 2
  }
  quiet <- function(V) {
    withCallingHandlers(svyrcs:::warn_if_rank_deficient(V, "x"),
                        warning = function(w) invokeRestart("muffleWarning"))
  }
  err <- function(V, beta) {
    e <- eigen(V, symmetric = TRUE)
    z <- crossprod(e$vectors, beta)
    qf <- as.numeric(t(beta) %*% solve(V) %*% beta)
    abs(qf - sum(z^2 / e$values)) / abs(sum(z^2 / e$values))
  }

  V7 <- make_cov(4, 1e7)
  V9 <- make_cov(4, 1e9)
  ## beta on the smallest eigenvector is the worst case: that is the direction the inverse amplifies
  worst <- function(V) eigen(V, symmetric = TRUE)$vectors[, 4]

  ## the guard has already given way by 1e9, while the arithmetic is still accurate
  expect_true(quiet(V9))
  expect_lt(err(V9, worst(V9)), 1e-6)

  ## and it stays quiet where the arithmetic is exact
  expect_false(quiet(make_cov(4, 1e2)))
  expect_lt(err(V7, worst(V7)), 1e-6)

  ## solve() itself holds far longer than either
  expect_false(inherits(try(solve(make_cov(4, 1e14)), silent = TRUE), "try-error"))
})
