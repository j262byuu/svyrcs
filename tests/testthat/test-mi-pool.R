mi_fits <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      kn <- rcs_knots(nhanes_bmi$bmi, 4)
      f <- stats::as.formula(do.call(substitute, list(
        str2lang("bmi ~ rcs(age, knots = KN) + sex + tchol + hdl"), list(KN = kn)
      )))
      cached <<- component_fits(mi_design(), f)
    }
    cached
  }
})

test_that("pooling reproduces mitools::MIcombine", {
  skip_if_not_installed("mitools")
  fits <- mi_fits()
  pooled <- pool_fits(fits, survey::degf(mi_design()$designs[[1]]))
  ref <- mitools::MIcombine(fits)

  expect_equal(unname(coef(pooled)), unname(coef(ref)), tolerance = 1e-10)
  expect_equal(unname(vcov(pooled)), unname(vcov(ref)), tolerance = 1e-10)
  expect_equal(names(coef(pooled)), names(coef(fits[[1]])))
})

test_that("identical imputations give zero between-imputation variance", {
  kn <- rcs_knots(nhanes_bmi$bmi, 4)
  f <- stats::as.formula(do.call(substitute, list(
    str2lang("bmi ~ rcs(age, knots = KN) + sex"), list(KN = kn)
  )))
  fits <- component_fits(mi_design_degenerate(), f)
  pooled <- pool_fits(fits, survey::degf(mi_design_degenerate()$designs[[1]]))

  expect_equal(max(abs(pooled$B)), 0)
  expect_equal(pooled$variance, pooled$Ubar)
  expect_equal(coef(pooled), coef(fits[[1]]))
})

test_that("the between-imputation variance is real when imputations differ", {
  pooled <- pool_fits(mi_fits(), 31)
  expect_gt(max(abs(pooled$B)), 0)
  ## total variance must exceed the within-imputation variance on the diagonal
  expect_true(all(diag(pooled$variance) >= diag(pooled$Ubar)))
})

test_that("Barnard-Rubin degrees of freedom behave", {
  nu <- 31

  ## exact at zero between-imputation variance
  expect_identical(barnard_rubin_df(u = 1, b = 0, m = 5, nu_com = nu), nu)
  expect_identical(barnard_rubin_df(u = c(1, 2), b = c(0, 0), m = 5, nu_com = nu), c(nu, nu))

  ## never above the complete-data df, and decreasing in b
  bs <- c(1e-8, 1e-4, 1e-2, 0.1, 1)
  got <- barnard_rubin_df(u = rep(1, length(bs)), b = bs, m = 5, nu_com = nu)
  expect_true(all(got <= nu))
  expect_true(all(diff(got) < 0))

  ## and strictly below mitools' uncorrected rule, which is what makes the correction necessary
  r <- (1 + 1 / 5) * bs / 1
  mitools_df <- (5 - 1) * (1 + 1 / r)^2
  expect_true(all(got < mitools_df))
})

test_that("the pooled fit answers the model generics", {
  fits <- mi_fits()
  pooled <- pool_fits(fits, 31)

  expect_s3_class(pooled, "svyrcs_mifit")
  expect_true(inherits(pooled, class(fits[[1]])[1]))
  expect_equal(terms(pooled), terms(fits[[1]]))
  expect_equal(family(pooled)$family, family(fits[[1]])$family)
  expect_equal(nobs(pooled), nobs(fits[[1]]))
  expect_equal(nrow(model.frame(pooled)), nrow(model.frame(fits[[1]])))

  ## and the pieces the single-design machinery needs
  expect_equal(find_rcs_term(pooled)$var, "age")
  expect_equal(effect_measure(pooled)$measure, "Difference")
})

test_that("predicting from the pooled fit is refused", {
  pooled <- pool_fits(mi_fits(), 31)
  expect_error(predict(pooled, newdata = nhanes_bmi[1:3, ]), class = "svyrcs_error")
  expect_error(predict(pooled), "component fit")
})

test_that("degenerate and malformed imputation sets are rejected", {
  expect_error(pool_fits(mi_fits()[1], 31), class = "svyrcs_error")

  single <- survey::svydesign(id = ~psu, strata = ~strata, weights = ~weight, nest = TRUE,
                              data = mitools::imputationList(list(nhanes_bmi)))
  expect_error(mi_components(single), class = "svyrcs_error")
  expect_error(mi_components(single), "at least 2")
})

test_that("component designs must share their degrees of freedom", {
  designs <- mi_design()
  expect_equal(shared_degf(designs$designs), survey::degf(designs$designs[[1]]))

  mixed <- designs$designs
  mixed[[2]] <- subset(mixed[[2]], cycle == "2005-2006")
  expect_error(shared_degf(mixed), class = "svyrcs_error")
  expect_error(shared_degf(mixed), "different degrees of freedom")
})

test_that("an imputed design is recognised", {
  expect_true(is_mi_design(mi_design()))
  expect_false(is_mi_design(nhanes_design()))
  expect_false(is_mi_design(nhanes_bmi))
})

test_that("the fraction of missing information is between 0 and 1", {
  expect_equal(fraction_missing(u = 1, b = 0, m = 5), 0)
  fmi <- fraction_missing(u = 1, b = c(0.01, 0.1, 1), m = 5)
  expect_true(all(fmi > 0 & fmi < 1))
  expect_true(all(diff(fmi) > 0))
})
