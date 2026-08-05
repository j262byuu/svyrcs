test_that("the basis is numerically identical to rms::rcs for 3 to 7 knots", {
  skip_if_not_installed("rms")
  set.seed(42)
  for (x in list(rnorm(300, 25, 6), rexp(300, 0.2), runif(300, 0, 1))) {
    for (nk in 3:7) {
      kn <- rcs_knots(x, nk)
      mine <- bare(rcspline(x, kn))
      theirs <- rms::rcs(x, kn)
      attributes(theirs) <- list(dim = dim(theirs))
      expect_equal(mine, theirs, tolerance = 1e-10,
                   info = paste("nk =", nk))
    }
  }
})

test_that("a knot count and the equivalent explicit knots give the same basis", {
  set.seed(1)
  x <- rnorm(200, 25, 5)
  for (nk in 3:7) {
    expect_equal(bare(rcspline(x, nk)), bare(rcspline(x, rcs_knots(x, nk))))
  }
})

test_that("knot counts use Harrell's recommended quantiles", {
  set.seed(1)
  x <- rnorm(500)
  expect_equal(rcs_knots(x, 4),
               unname(quantile(x, c(0.05, 0.35, 0.65, 0.95))))
  expect_equal(rcs_knots(x, 3),
               unname(quantile(x, c(0.10, 0.50, 0.90))))
})

test_that("the basis has k - 1 columns and carries its knots", {
  set.seed(2)
  x <- rnorm(100, 10, 2)
  for (nk in 3:7) {
    b <- rcspline(x, nk)
    expect_s3_class(b, "svyrcs_basis")
    expect_equal(ncol(b), nk - 1L)
    expect_equal(length(attr(b, "knots")), nk)
    expect_equal(attr(b, "nk"), nk)
    ## the first column is the exposure itself, so its coefficient is the linear term
    expect_equal(as.numeric(b[, 1]), x)
  }
})

test_that("row subsetting preserves the knots", {
  set.seed(3)
  b <- rcspline(rnorm(100, 10, 2), 4)
  b2 <- b[1:10, ]
  expect_s3_class(b2, "svyrcs_basis")
  expect_equal(attr(b2, "knots"), attr(b, "knots"))
  expect_equal(nrow(b2), 10L)
})

test_that("predicting on new data reuses the fitted knots", {
  set.seed(4)
  d <- data.frame(x = rnorm(300, 25, 6))
  d$y <- 0.1 * d$x + rnorm(300)
  fit <- lm(y ~ rcspline(x, 4), data = d)

  fitted_knots <- rcs_knots(d$x, 4)
  ## a grid whose own quantiles are nothing like the training data's
  nd <- data.frame(x = c(5, 10, 15))
  mf <- model.frame(delete.response(terms(fit)), nd)
  basis_new <- mf[["rcspline(x, 4)"]]

  expect_equal(bare(basis_new), bare(rcs_design_matrix(nd$x, fitted_knots)))
  ## and the whole point: three points cannot themselves support four data-driven knots
  expect_error(rcspline(nd$x, 4), class = "svyrcs_error")
})

test_that("invalid knot specifications are rejected with an svyrcs error", {
  set.seed(5)
  x <- rnorm(100)
  expect_error(rcspline(x, 2), class = "svyrcs_error")
  expect_error(rcspline(x, 8), class = "svyrcs_error")
  expect_error(rcspline(x, 4.5), class = "svyrcs_error")
  expect_error(rcspline(letters), class = "svyrcs_error")
  expect_error(rcspline(x, c(1, 2)), class = "svyrcs_error")
  expect_error(rcspline(x, c(1, 2, NA)), class = "svyrcs_error")
  expect_error(rcspline(rep(1, 10), 4), class = "svyrcs_error")
  expect_error(rcspline(rep(NA_real_, 10), 4), class = "svyrcs_error")
})

test_that("explicit knots are sorted and de-duplicated", {
  set.seed(6)
  x <- rnorm(200, 25, 5)
  expect_equal(rcs_knots(x, c(30, 20, 25)), c(20, 25, 30))
  expect_equal(rcs_knots(x, c(20, 25, 25, 30)), c(20, 25, 30))
})

test_that("coefficient names are the term label with primes", {
  set.seed(7)
  d <- data.frame(x = rnorm(200, 25, 5), z = rnorm(200))
  d$y <- rnorm(200)
  fit <- lm(y ~ rcspline(x, 4) + z, data = d)
  expect_equal(names(coef(fit)),
               c("(Intercept)", "rcspline(x, 4)", "rcspline(x, 4)'", "rcspline(x, 4)''", "z"))
})
