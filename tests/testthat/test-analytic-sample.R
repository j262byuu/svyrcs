## Review item 1: knots, reference, range and the extrapolation check must come from the rows the
## model is actually fitted on, not from the design as supplied. Review items 7a and 7b: formula
## forms that were rejected or mis-diagnosed.

analytic_of <- function(design, vars) design[stats::complete.cases(nhanes_bmi[vars]), ]

test_that("missingness related to the exposure moves the knots, and the fit follows the analytic sample", {
  ## Drop most of the high-BMI rows, so the two samples genuinely disagree.
  design <- nhanes_design()
  set.seed(1)
  hi <- which(nhanes_bmi$bmi > 32)
  design$variables$y <- nhanes_bmi$tchol
  design$variables$y[hi[seq_len(round(0.9 * length(hi)))]] <- NA

  fit <- svyrcs(y ~ rcspline(bmi, 4) + age, design = design, n = 20)
  keep <- stats::complete.cases(design$variables[c("y", "bmi", "age")])
  analytic <- design[keep, ]

  expect_equal(fit$n, sum(keep))
  expect_equal(fit$knots,
               signif(as.numeric(coef(survey::svyquantile(~bmi, analytic, c(0.05, 0.35, 0.65, 0.95),
                                                          na.rm = TRUE, ci = FALSE))), 7))
  expect_equal(fit$ref$value,
               as.numeric(coef(survey::svyquantile(~bmi, analytic, 0.5, na.rm = TRUE, ci = FALSE))))

  ## and the full-design values, which the old implementation used, are far away
  full_kn <- as.numeric(coef(survey::svyquantile(~bmi, design, c(0.05, 0.35, 0.65, 0.95),
                                                 na.rm = TRUE, ci = FALSE)))
  expect_gt(abs(fit$knots[4] - full_kn[4]), 5)
})

test_that("the plotting range and the extrapolation check use the analytic sample", {
  design <- nhanes_design()
  set.seed(2)
  hi <- which(nhanes_bmi$bmi > 32)
  design$variables$y <- nhanes_bmi$tchol
  design$variables$y[hi[seq_len(round(0.9 * length(hi)))]] <- NA
  fit <- svyrcs(y ~ rcspline(bmi, 4) + age, design = design, n = 20)
  analytic <- design[stats::complete.cases(design$variables[c("y", "bmi", "age")]), ]

  expect_equal(range(fit$curve$x),
               as.numeric(coef(survey::svyquantile(~bmi, analytic, c(0.01, 0.99),
                                                   na.rm = TRUE, ci = FALSE))))

  ## a point inside the full design's range but outside the analytic sample's must now warn
  outside <- max(analytic$variables$bmi[!is.na(analytic$variables$y)], na.rm = TRUE) + 20
  expect_lt(outside, max(nhanes_bmi$bmi))
  w <- capture_warnings(predict(fit, x = outside))
  expect_true(any(grepl("outside the observed range", w)))
})

test_that("with no missing data nothing changes", {
  ## bmi, age and sex are complete, so the analytic sample is the whole design and the numbers must
  ## be exactly what they were.
  design <- nhanes_design()
  fit <- svyrcs(bmi ~ rcspline(age, 4) + sex, design = design, n = 20)
  expect_equal(fit$n, nrow(nhanes_bmi))
  expect_equal(fit$n_dropped, 0L)
  expect_equal(fit$knots,
               signif(as.numeric(coef(survey::svyquantile(~age, design, c(0.05, 0.35, 0.65, 0.95),
                                                          na.rm = TRUE, ci = FALSE))), 7))
})

test_that("the mask covers every variable in the formula, not only the exposure", {
  design <- nhanes_design()
  design$variables$cov <- nhanes_bmi$hdl        # missing where hdl is
  design$variables$cov[1:800] <- NA             # and some more
  fit <- svyrcs(bmi ~ rcspline(age, 4) + cov, design = design, n = 5)
  expect_equal(fit$n, sum(stats::complete.cases(design$variables[c("bmi", "age", "cov")])))
  expect_gt(fit$n_dropped, 800L)
})

test_that("under imputation every imputation is fitted on the same rows", {
  ii <- lapply(1:4, function(i) {
    set.seed(i)
    o <- nhanes_bmi
    m <- is.na(o$tchol)
    o$tchol[m] <- sample(o$tchol[!m], sum(m), TRUE)
    ## outcome missingness depends on the imputed value, so the per-imputation masks differ
    o$y <- ifelse(o$tchol > 260, NA, o$bmi)
    o
  })
  des <- survey::svydesign(id = ~psu, strata = ~strata, weights = ~weight, nest = TRUE,
                           data = mitools::imputationList(ii))
  fit <- svyrcs(y ~ rcspline(age, 4) + tchol, design = des, n = 10)

  masks <- lapply(ii, function(o) stats::complete.cases(o[c("y", "age", "tchol")]))
  shared <- Reduce(`&`, masks)
  expect_equal(fit$n, sum(shared))
  expect_lt(sum(shared), max(vapply(masks, sum, integer(1L))))
  expect_equal(fit$imputations$shared_mask_cost,
               max(vapply(masks, sum, integer(1L))) - sum(shared))

  ## every component fit really did use the same number of rows
  ns <- vapply(fit$model$fits, function(f) as.integer(stats::nobs(f)), integer(1L))
  expect_true(all(ns == sum(shared)))
  expect_output(print(fit), "shares one sample")
})

test_that("a namespaced rcspline() call is recognised", {
  design <- nhanes_design()
  a <- svyrcs(tchol ~ rcspline(bmi, 4) + age, design = design, n = 10)
  b <- svyrcs(tchol ~ svyrcs::rcspline(bmi, 4) + age, design = design, n = 10)
  expect_equal(b$var, "bmi")
  expect_equal(b$knots, a$knots)
  expect_equal(b$curve$estimate, a$curve$estimate)

  ## and it works inside an interaction too
  d <- svyrcs(tchol ~ svyrcs::rcspline(bmi, 4) * sex + age, design = design, n = 10)
  expect_equal(d$groups$var, "sex")
})

test_that("a transformed exposure is refused with an explanation", {
  design <- nhanes_design()
  err <- tryCatch(svyrcs(tchol ~ rcspline(log(bmi), 4) + age, design = design),
                  svyrcs_error = function(e) conditionMessage(e))
  expect_match(err, "bare variable name")
  expect_match(err, "log\\(bmi\\)")
  expect_match(err, "update\\(design")

  ## the documented workaround works
  design$variables$log_bmi <- log(nhanes_bmi$bmi)
  expect_s3_class(svyrcs(tchol ~ rcspline(log_bmi, 4) + age, design = design, n = 5), "svyrcs")
})
