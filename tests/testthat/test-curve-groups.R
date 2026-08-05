test_that("a grouped curve has one block per level, in level order", {
  fit <- int_model("rcspline(bmi, knots = KN) * race + age")
  cv <- svyrcs_curve(fit, "bmi", n = 20)

  expect_true("group" %in% names(cv))
  expect_s3_class(cv$group, "factor")
  expect_equal(levels(cv$group), levels(nhanes_bmi$race))
  expect_equal(nrow(cv), 20L * nlevels(nhanes_bmi$race))
  expect_equal(as.integer(table(cv$group)), rep(20L, nlevels(nhanes_bmi$race)))
  expect_true(all(is.finite(cv$estimate)))
  expect_equal(attr(cv, "modifier")$var, "race")
})

test_that("every group's estimate at the reference is exactly the null with zero SE", {
  fit <- int_model("rcspline(bmi, knots = KN) * sex + age", cox = TRUE,
                   response = "survival::Surv(time, event)")
  cv <- svyrcs_curve(fit, "bmi", ref = 27)
  at_ref <- svyrcs_curve(fit, "bmi", ref = 27, at = 27)
  expect_identical(unique(at_ref$estimate), 1)
  expect_identical(unique(at_ref$se), 0)
  expect_equal(nrow(at_ref), 2L)
  expect_equal(attr(cv, "ref"), 27)
})

test_that("per-group contrasts reproduce predict(type = 'lp') within each subgroup", {
  skip_if_not_installed("survival")
  fit <- int_model("rcspline(bmi, knots = KN) * sex + age", cox = TRUE,
                   response = "survival::Surv(time, event)")
  design <- nhanes_design()
  grid <- seq(18, 45, length.out = 25)
  x0 <- 27
  cv <- svyrcs_curve(fit, "bmi", ref = x0, at = grid)

  for (lvl in levels(nhanes_bmi$sex)) {
    row <- which(design$variables$sex == lvl)[1L]
    nd <- design$variables[rep(row, length(grid)), ]
    nd$bmi <- grid
    lp <- as.numeric(predict(fit, newdata = nd, type = "lp"))
    nd0 <- nd[1L, ]
    nd0$bmi <- x0
    lp0 <- as.numeric(predict(fit, newdata = nd0, type = "lp"))

    got <- cv[cv$group == lvl, ]
    expect_equal(log(got$estimate), lp - lp0, tolerance = 1e-8, info = lvl)
  }
})

test_that("both formula orders give identical curves", {
  a <- svyrcs_curve(int_model("rcspline(bmi, knots = KN) * sex + age"), "bmi", ref = 27, n = 15)
  b <- svyrcs_curve(int_model("sex * rcspline(bmi, knots = KN) + age"), "bmi", ref = 27, n = 15)
  expect_equal(as.data.frame(a)[c("x", "estimate", "conf.low", "conf.high", "se")],
               as.data.frame(b)[c("x", "estimate", "conf.low", "conf.high", "se")])
  expect_equal(as.character(a$group), as.character(b$group))
})

test_that("the reference group's curve equals the same model fitted without interaction terms", {
  ## The reference level uses main effects only, so its estimates must be exactly what the
  ## selection matrix extracts from the full coefficient vector.
  fit <- int_model("rcspline(bmi, knots = KN) * sex + age")
  term <- find_rcs_term(fit)
  m <- find_modifier(fit, term)
  b <- coef(fit)
  V <- vcov(fit)
  kn <- term$knots

  x <- c(20, 30, 40)
  dB <- sweep(rcs_design_matrix(x, kn), 2, as.numeric(rcs_design_matrix(27, kn)))
  manual <- as.numeric(dB %*% b[m$spline_cols])

  got <- svyrcs_curve(fit, "bmi", ref = 27, at = x, group = "Male")
  expect_equal(got$estimate, manual, tolerance = 1e-10)
})

test_that("a non-reference group's variance includes the main-interaction covariance", {
  fit <- int_model("rcspline(bmi, knots = KN) * sex + age")
  term <- find_rcs_term(fit)
  m <- find_modifier(fit, term)
  b <- coef(fit)
  V <- vcov(fit)
  kn <- term$knots
  x <- 40
  dB <- sweep(rcs_design_matrix(x, kn), 2, as.numeric(rcs_design_matrix(27, kn)))

  got <- svyrcs_curve(fit, "bmi", ref = 27, at = x, group = "Female")

  main <- m$spline_cols
  inter <- m$columns$Female
  L <- group_selection(names(b), term, m, "Female")
  correct <- as.numeric(rowSums(((dB %*% L) %*% V) * (dB %*% L)))
  expect_equal(got$se^2, correct, tolerance = 1e-12)

  ## the naive "add the two variances" answer is different, which is the point of the matrix form
  naive <- as.numeric(rowSums((dB %*% V[main, main]) * dB)) +
    as.numeric(rowSums((dB %*% V[inter, inter]) * dB))
  expect_false(isTRUE(all.equal(correct, naive)))
})

test_that("group selects a subset of levels and rejects unknown ones", {
  fit <- int_model("rcspline(bmi, knots = KN) * race + age")
  one <- svyrcs_curve(fit, "bmi", n = 10, group = "Other")
  expect_equal(nrow(one), 10L)
  expect_equal(as.character(unique(one$group)), "Other")

  two <- svyrcs_curve(fit, "bmi", n = 10, group = c("Hispanic", "Other"))
  expect_equal(nrow(two), 20L)

  expect_error(svyrcs_curve(fit, "bmi", group = "Martian"), class = "svyrcs_error")
  expect_error(svyrcs_curve(bare_glm(), "bmi", group = "Male"), class = "svyrcs_error")
})

test_that("a curve-derived reference is located on the reference level and says so", {
  skip_if_not_installed("survival")
  fit <- int_model("rcspline(bmi, knots = KN) * sex + age", cox = TRUE,
                   response = "survival::Surv(time, event)")
  cv <- svyrcs_curve(fit, "bmi", ref = "min", n = 300)
  expect_match(attr(cv, "ref_method"), "minimum-risk point \\(Male\\)")

  ## anchored at its own minimum, the reference level's curve cannot dip below the null
  male <- cv[cv$group == "Male", ]
  expect_gte(min(male$estimate), 1 - 1e-8)
})

test_that("all groups share one reference value", {
  fit <- int_model("rcspline(bmi, knots = KN) * race + age")
  cv <- svyrcs_curve(fit, "bmi", n = 50)
  ref <- attr(cv, "ref")
  at_ref <- svyrcs_curve(fit, "bmi", at = ref)
  expect_equal(at_ref$estimate, rep(0, nlevels(nhanes_bmi$race)))
  expect_equal(at_ref$se, rep(0, nlevels(nhanes_bmi$race)))
})
