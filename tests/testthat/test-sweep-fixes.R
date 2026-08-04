## Regression tests for the adversarial sweep of 2026-08-03. Each one uses the input that exposed
## the problem, so a reintroduction fails here rather than in someone's analysis.

test_that("an ordered factor works as a modifier", {
  ## Ordered factors use polynomial contrasts, so their columns are `agegrp.L` and `agegrp.Q` and a
  ## level name never appears in a coefficient name. Constructing columns from level names could
  ## therefore never match.
  design <- nhanes_design()
  design$variables$agegrp <- ordered(cut(nhanes_bmi$age, c(0, 40, 60, 100)))

  fit <- svyrcs(tchol ~ rcs(bmi, 4) * agegrp + sex, design = design, n = 20)
  expect_equal(fit$groups$var, "agegrp")
  expect_equal(length(fit$groups$levels), 3L)
  expect_equal(nrow(fit$curve), 60L)
  expect_true(all(is.finite(fit$curve$estimate)))
  expect_true(is.finite(fit$tests$interaction$p_F))
  expect_equal(fit$tests$interaction$df1, 3L * 2L)
})

test_that("an ordered and an unordered modifier give the same fitted curves", {
  ## The contrast coding changes the parameterisation, not the fit, so the per-group curves must
  ## agree even though the coefficients do not.
  design <- nhanes_design()
  cuts <- cut(nhanes_bmi$age, c(0, 40, 60, 100))
  design$variables$ag_ord <- ordered(cuts)
  design$variables$ag_unord <- factor(as.character(cuts))

  a <- svyrcs(tchol ~ rcs(bmi, 4) * ag_ord + sex, design = design, n = 20)
  b <- svyrcs(tchol ~ rcs(bmi, 4) * ag_unord + sex, design = design, n = 20)

  expect_equal(a$curve$estimate, b$curve$estimate, tolerance = 1e-8)
  expect_equal(a$curve$se, b$curve$se, tolerance = 1e-8)
  expect_equal(a$tests$interaction$F, b$tests$interaction$F, tolerance = 1e-8)
})

test_that("the reference level of a treatment-contrast modifier still uses main effects only", {
  ## Guards the generalisation: with treatment contrasts the contrast row is 0 for the reference
  ## level and a single 1 elsewhere, so the selection matrix must be unchanged.
  fit <- svyrcs(tchol ~ rcs(bmi, 4) * sex + age, design = nhanes_design(), n = 10)
  term <- find_rcs_term(fit$model)
  m <- find_modifier(fit$model, term)
  b <- coef(fit$model)

  L_ref <- group_selection(names(b), term, m, m$ref_level)
  expect_equal(as.numeric(L_ref %*% b), unname(b[m$spline_cols]))
  expect_equal(sort(unique(as.numeric(L_ref))), c(0, 1))

  L_oth <- group_selection(names(b), term, m, "Female")
  expect_equal(as.numeric(L_oth %*% b), unname(b[m$spline_cols] + b[m$columns$Female]))
})

test_that("a rank-deficient design warns on both entry points instead of erroring on one", {
  ## Previously svyrcs() stopped while svyrcs_curve() returned a band for the same data.
  design <- nhanes_design()
  tiny <- subset(design, strata %in% sort(unique(nhanes_bmi$strata))[1:2])
  expect_lte(survey::degf(tiny), 2)

  ## Several warnings fire here -- rank deficiency from the curve, unavailable tests from the Wald
  ## block -- so assert on the whole set rather than the first.
  w <- capture_warnings(fit <- svyrcs(tchol ~ rcs(bmi, 4) + age, design = tiny, n = 5))
  expect_true(any(grepl("cannot identify this model", w)))
  expect_true(any(grepl("singular", w)))
  expect_true(is.na(fit$tests$overall$p_F))
  expect_true(all(is.finite(fit$curve$estimate)))

  kn <- fit$knots
  f <- stats::as.formula(do.call(substitute, list(
    str2lang("tchol ~ rcs(bmi, knots = KN) + age"), list(KN = kn)
  )))
  m <- survey::svyglm(f, design = tiny)
  expect_warning(svyrcs_curve(m, "bmi", n = 5), "rank")
})

test_that("a full-rank design warns about nothing", {
  expect_silent(svyrcs(tchol ~ rcs(bmi, 4) + age, design = nhanes_design(), n = 5))
})

test_that("evaluating outside the observed exposure range warns", {
  fit <- svyrcs(tchol ~ rcs(bmi, 4) + age, design = nhanes_design(), n = 5)
  obs <- range(nhanes_bmi$bmi)

  w <- capture_warnings(predict(fit, x = c(-50, 500)))
  expect_true(any(grepl("outside the observed range", w)))
  expect_true(any(grepl("extrapolation", w)))
  w2 <- capture_warnings(svyrcs(tchol ~ rcs(bmi, 4) + age, design = nhanes_design(),
                                range = c(200, 300), n = 5))
  expect_true(any(grepl("outside the observed range", w2)))

  ## inside the range, nothing is said
  expect_silent(predict(fit, x = c(obs[1], mean(obs), obs[2])))

  ## a reference outside the range is still refused outright, because it would make every estimate
  ## an extrapolation rather than just the point asked for
  expect_error(predict(fit, x = 30, ref = 500), class = "svyrcs_error")
})

test_that("colliding weighted quantiles say what actually went wrong", {
  design <- nhanes_design()
  design$variables$few <- as.numeric(cut(nhanes_bmi$bmi, 5, labels = FALSE))
  err <- tryCatch(svyrcs(tchol ~ rcs(few, 4) + age, design = design),
                  svyrcs_error = function(e) conditionMessage(e))
  expect_match(err, "could not place 4 distinct knots")
  expect_match(err, "quantiles coincide")
  expect_match(err, "few distinct values")
})

test_that("imputations with different coefficients name the culprit", {
  ii <- lapply(1:3, function(i) {
    set.seed(i); o <- nhanes_bmi
    mm <- is.na(o$tchol); o$tchol[mm] <- sample(o$tchol[!mm], sum(mm), TRUE)
    o
  })
  ii[[2]]$race <- droplevels(factor(ifelse(ii[[2]]$race == "Other", "Hispanic",
                                           as.character(ii[[2]]$race))))
  des <- survey::svydesign(id = ~psu, strata = ~strata, weights = ~weight, nest = TRUE,
                           data = mitools::imputationList(ii))
  err <- tryCatch(svyrcs(bmi ~ rcs(age, 4) + race, design = des),
                  svyrcs_error = function(e) conditionMessage(e))
  expect_match(err, "different coefficients")
  expect_match(err, "imputation 2")
})

test_that("a boundary extremum is labelled as such", {
  ## On a monotone curve the minimum-risk point is just the end of the plotting range, which is a
  ## property of `range` rather than a turning point in the data.
  fit <- svyrcs(tchol ~ rcs(age, 4) + sex, design = nhanes_design(), ref = "min", n = 100)
  if (fit$ref$value %in% range(fit$curve$x)) {
    expect_match(fit$ref$method, "at the range boundary")
  } else {
    expect_false(grepl("range boundary", fit$ref$method))
  }
})

test_that("a log link on a gaussian family is a ratio, not a risk ratio", {
  fit <- svyrcs(tchol ~ rcs(bmi, 4) + age, design = nhanes_design(),
                family = gaussian(link = "log"), n = 5)
  expect_equal(fit$measure, "Ratio")
  expect_equal(fit$null, 1)

  ## A count model without person-time compares expected counts, not rates.
  fit2 <- svyrcs(statin_use ~ rcs(bmi, 4) + age, design = nhanes_design(),
                 family = quasipoisson(), n = 5)
  expect_equal(fit2$measure, "Mean ratio")

  ## with an offset it really is a rate ratio
  design <- nhanes_design()
  design$variables$py <- pmax(nhanes_bmi$time, 0.1)
  fit3 <- svyrcs(event ~ rcs(bmi, 4) + age + offset(log(py)), design = design,
                 family = quasipoisson(), n = 5)
  expect_equal(fit3$measure, "Rate ratio")
})

test_that("non-finite estimates are flagged rather than returned silently", {
  design <- nhanes_design()
  design$variables$sep <- as.integer(nhanes_bmi$bmi > 30)
  ## Perfect separation also makes the fit rank deficient, so more than one warning is raised;
  ## what matters is that the non-finite estimates are among the things reported.
  w <- capture_warnings(
    svyrcs(sep ~ rcs(bmi, 4) + age, design = design, family = quasibinomial(), n = 5)
  )
  expect_true(any(grepl("non-finite", w)))
})
