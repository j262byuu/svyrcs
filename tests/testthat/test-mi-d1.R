## D1 (Li, Raghunathan and Rubin 1991) with Reiter (2007) small-sample degrees of freedom.
## mitml is the reference implementation; where both are defined they must agree exactly.

d1_inputs <- function(fit, cols) {
  list(
    Qbar = coef(fit$model)[cols],
    Ubar = fit$model$Ubar[cols, cols, drop = FALSE],
    B = fit$model$B[cols, cols, drop = FALSE],
    m = fit$model$m
  )
}

oracle_inputs <- function(fit, cols) {
  m <- fit$model$m
  k <- length(cols)
  Q <- vapply(fit$model$fits, function(x) coef(x)[cols], numeric(k))
  U <- array(vapply(fit$model$fits, function(x) vcov(x)[cols, cols], numeric(k * k)),
             dim = c(k, k, m))
  list(Qhat = Q, Uhat = U)
}

test_that("d1_test agrees with mitml across m and block size", {
  skip_if_not_installed("mitml")
  oracle <- getFromNamespace(".D1", "mitml")

  for (m in c(5L, 10L)) {
    des <- survey::svydesign(
      id = ~psu, strata = ~strata, weights = ~weight, nest = TRUE,
      data = mitools::imputationList(deterministic_imputations(m, seed = 3))
    )
    fit <- svyrcs(bmi ~ rcspline(age, 4) + tchol + hdl, design = des, n = 5)
    term <- find_rcs_term(fit$model)
    cols <- spline_colnames(term$label, term$nk)

    for (block in list(cols, cols[-1L])) {
      k <- length(block)
      skip_if(k * (m - 1) <= 4)
      inp <- d1_inputs(fit, block)
      mine <- d1_test(inp$Qbar, inp$Ubar, inp$B, inp$m, fit$degf)
      oi <- oracle_inputs(fit, block)
      ref <- oracle(Qhat = oi$Qhat, Uhat = oi$Uhat, df.com = fit$degf)

      expect_equal(mine$F, as.numeric(ref$F), tolerance = 1e-10,
                   info = paste("m", m, "k", k))
      expect_equal(mine$v, ref$v, tolerance = 1e-10, info = paste("m", m, "k", k))
      expect_equal(mine$r, ref$r, tolerance = 1e-12, info = paste("m", m, "k", k))
    }
  }
})

test_that("the reported joint tests are the D1 numbers", {
  des <- survey::svydesign(
    id = ~psu, strata = ~strata, weights = ~weight, nest = TRUE,
    data = mitools::imputationList(deterministic_imputations(5, seed = 4))
  )
  fit <- svyrcs(bmi ~ rcspline(age, 4) + tchol, design = des, n = 5)
  term <- find_rcs_term(fit$model)
  cols <- spline_colnames(term$label, term$nk)

  inp <- d1_inputs(fit, cols)
  d1 <- d1_test(inp$Qbar, inp$Ubar, inp$B, inp$m, fit$degf)
  expect_equal(fit$tests$overall$F, d1$F, tolerance = 1e-12)
  expect_equal(fit$tests$overall$df2, d1$v, tolerance = 1e-12)
  expect_equal(fit$tests$overall$fmi, d1$r / (1 + d1$r), tolerance = 1e-12)
  expect_equal(fit$tests$overall$p_F, pf(d1$F, length(cols), d1$v, lower.tail = FALSE),
               tolerance = 1e-12)
})

test_that("the interaction block goes through D1 too", {
  skip_if_not_installed("mitml")
  oracle <- getFromNamespace(".D1", "mitml")
  des <- survey::svydesign(
    id = ~psu, strata = ~strata, weights = ~weight, nest = TRUE,
    data = mitools::imputationList(deterministic_imputations(5, seed = 5))
  )
  fit <- svyrcs(bmi ~ rcspline(age, 4) * sex + tchol, design = des, n = 5)
  m <- find_modifier(fit$model, find_rcs_term(fit$model))
  icols <- unlist(m$columns, use.names = FALSE)

  oi <- oracle_inputs(fit, icols)
  ref <- oracle(Qhat = oi$Qhat, Uhat = oi$Uhat, df.com = fit$degf)
  expect_equal(fit$tests$interaction$F, as.numeric(ref$F), tolerance = 1e-10)
  expect_equal(fit$tests$interaction$df2, ref$v, tolerance = 1e-10)
})

test_that("df reaches the imputation arithmetic", {
  ## In 0.5.1 `df` affected only the single-design path, so this silently did nothing.
  des <- survey::svydesign(
    id = ~psu, strata = ~strata, weights = ~weight, nest = TRUE,
    data = mitools::imputationList(deterministic_imputations(5, seed = 6))
  )
  base <- svyrcs(bmi ~ rcspline(age, 4) + tchol, design = des, n = 5)
  small <- svyrcs(bmi ~ rcspline(age, 4) + tchol, design = des, n = 5, df = 10)
  big <- svyrcs(bmi ~ rcspline(age, 4) + tchol, design = des, n = 5, df = Inf)

  expect_lt(small$tests$overall$df2, base$tests$overall$df2)
  expect_gt(big$tests$overall$df2, base$tests$overall$df2)
  expect_lt(max(small$curve$df), max(base$curve$df))
  expect_gt(min(big$curve$df), max(base$curve$df))

  ## the statistic itself does not depend on the reference distribution
  expect_equal(small$tests$overall$F, base$tests$overall$F)
  expect_equal(big$tests$overall$F, base$tests$overall$F)
})

test_that("a singular within-imputation covariance warns rather than stopping", {
  Ubar <- matrix(c(1, 1, 1, 1), 2, 2)
  expect_null(d1_test(c(1, 2), Ubar, matrix(0, 2, 2), m = 5, dfcom = 31))

  expect_warning(out <- wald_block(c(1, 2), Ubar, 31, "overall association",
                                   mi = list(Ubar = Ubar, B = matrix(0, 2, 2), m = 5, degf = 31)),
                 "singular")
  expect_true(is.na(out$p_F))
})
