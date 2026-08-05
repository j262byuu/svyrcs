## Validation, labelling and documented-contract regressions from the Codex review.

test_that("default knot placement matches Hmisc, including the small-sample rule", {
  skip_if_not_installed("Hmisc")
  ## Below 100 observations Harrell moves the outer knots to the fifth smallest and fifth largest
  ## values. The package previously used plain quantiles throughout, so the documented equivalence
  ## with rms held for the basis but not for the default placement.
  set.seed(9)
  for (n in c(10, 20, 50, 99, 100, 150, 300)) {
    for (nk in 3:7) {
      for (x in list(as.numeric(seq_len(n)), rnorm(n), rgamma(n, 3, 0.5))) {
        expect_equal(rcs_knots(x, nk),
                     as.numeric(Hmisc::rcspline.eval(x, nk = nk, knots.only = TRUE)),
                     info = paste("n", n, "nk", nk))
      }
    }
  }
})

test_that("the small-sample rule only applies below 100 observations", {
  x99 <- as.numeric(seq_len(99))
  x100 <- as.numeric(seq_len(100))
  probs <- harrell_knot_probs(4)

  ## at n = 99 the outer knots are the 5th smallest and 5th largest
  expect_equal(rcs_knots(x99, 4)[c(1, 4)], c(5, 95))
  ## at n = 100 they are the plain quantiles again
  expect_equal(rcs_knots(x100, 4), unname(quantile(x100, probs)))
})

test_that("knot probabilities are computed, not tabulated", {
  ## The published table rounds to four decimals, which put 7-knot placement about 0.004 away from
  ## Hmisc.
  for (nk in 3:7) {
    outer <- if (nk > 6) 0.025 else if (nk > 3) 0.05 else 0.1
    expect_equal(harrell_knot_probs(nk), seq(outer, 1 - outer, length.out = nk))
  }
  expect_error(harrell_knot_probs(2), class = "svyrcs_error")
  expect_error(harrell_knot_probs(8), class = "svyrcs_error")
})

test_that("non-finite knots are refused", {
  expect_error(rcspline(1:20, c(1, 2, Inf)), class = "svyrcs_error")
  expect_error(rcspline(1:20, c(1, NaN, 3)), class = "svyrcs_error")
  expect_error(rcspline(1:20, c(-Inf, 2, 3)), class = "svyrcs_error")
  expect_error(rcspline(1:20, Inf), class = "svyrcs_error")
  expect_error(rcspline(1:20, NA), class = "svyrcs_error")
  expect_error(rcspline(1:20, NaN), class = "svyrcs_error")

  ## and the message names finiteness rather than blaming something else
  expect_error(rcspline(1:20, c(1, 2, Inf)), "finite")
})

test_that("range and at are validated", {
  des <- nhanes_design()
  fit <- svyrcs(tchol ~ rcspline(bmi, 4) + age, design = des, n = 5)

  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = des, range = 30), class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = des, range = c(NA, 40)),
               class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = des, range = c(20, Inf)),
               class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = des, range = c(30, 30)),
               class = "svyrcs_error")

  expect_error(predict(fit, x = c(30, NA)), class = "svyrcs_error")
  expect_error(predict(fit, x = Inf), class = "svyrcs_error")
  expect_error(predict(fit, x = c(30, NaN)), class = "svyrcs_error")

  ## a reversed range is accepted and sorted, which is harmless
  rev <- svyrcs(tchol ~ rcspline(bmi, 4), design = des, range = c(40, 20), n = 5)
  expect_equal(range(rev$curve$x), c(20, 40))
})

test_that("level and ref_prob validation is unchanged", {
  ## The review claimed `level = NA` raised an untyped error. It does not, and this pins that.
  des <- nhanes_design()
  fit <- svyrcs(tchol ~ rcspline(bmi, 4) + age, design = des, n = 5)
  expect_error(predict(fit, x = 30, level = NA), class = "svyrcs_error")
  expect_error(predict(fit, x = 30, level = 1.5), class = "svyrcs_error")
  expect_error(svyrcs(tchol ~ rcspline(bmi, 4), design = des, ref = "quantile", ref_prob = NA),
               class = "svyrcs_error")
})

test_that("a log link is labelled by what the family estimates", {
  des <- nhanes_design()
  des$variables$py <- pmax(nhanes_bmi$time, 0.1)

  ## binomial: the fitted mean is a probability, so the contrast is a risk ratio however it is offset
  binom <- suppressWarnings(svyrcs(high_chol ~ rcspline(bmi, 4) + age, design = des,
                                   family = quasibinomial(link = "log"), n = 5))
  expect_equal(binom$measure, "RR")

  ## count family: person-time makes it a rate, its absence makes it a ratio of expected counts
  expect_equal(svyrcs(event ~ rcspline(bmi, 4) + age + offset(log(py)), design = des,
                      family = quasipoisson(), n = 5)$measure, "Rate ratio")
  expect_equal(svyrcs(statin_use ~ rcspline(bmi, 4) + age, design = des,
                      family = quasipoisson(), n = 5)$measure, "Mean ratio")

  ## anything else is just a ratio
  expect_equal(svyrcs(tchol ~ rcspline(bmi, 4) + age, design = des,
                      family = gaussian(link = "log"), n = 5)$measure, "Ratio")
})

test_that("every supported family produces a documented measure", {
  ## The @return block listed values the code had stopped producing. This asserts the contract so it
  ## cannot drift again without a failure.
  documented <- c("HR", "OR", "RR", "Rate ratio", "Mean ratio", "Ratio", "Difference")
  des <- nhanes_design()
  des$variables$py <- pmax(nhanes_bmi$time, 0.1)
  skip_if_not_installed("survival")

  fits <- list(
    svyrcs(survival::Surv(time, event) ~ rcspline(bmi, 4) + age, design = des, n = 5),
    svyrcs(high_chol ~ rcspline(bmi, 4) + age, design = des, family = quasibinomial(), n = 5),
    suppressWarnings(svyrcs(high_chol ~ rcspline(bmi, 4) + age, design = des,
                            family = quasibinomial(link = "log"), n = 5)),
    svyrcs(event ~ rcspline(bmi, 4) + age + offset(log(py)), design = des,
           family = quasipoisson(), n = 5),
    svyrcs(statin_use ~ rcspline(bmi, 4) + age, design = des, family = quasipoisson(), n = 5),
    svyrcs(tchol ~ rcspline(bmi, 4) + age, design = des, family = gaussian(link = "log"), n = 5),
    svyrcs(tchol ~ rcspline(bmi, 4) + age, design = des, n = 5)
  )
  measures <- vapply(fits, function(f) f$measure, character(1L))
  expect_true(all(measures %in% documented), info = paste(measures, collapse = ", "))
  expect_setequal(measures, documented)
})

test_that("the shipped documentation no longer contradicts the implementation", {
  ## Checked against the installed help rather than the README: the README is not installed, so
  ## reading it passed from a source checkout and failed under R CMD check. Aiming at the help also
  ## caught a second copy of the claim in ?svyrcs that grepping the README for its exact wording had
  ## missed.
  ## Under devtools::load_all() the installed help belongs to whatever version is in the library,
  ## which need not match the sources just loaded. R CMD check installs first, so there the database
  ## is the one being tested; anywhere else, skip rather than assert against the wrong package.
  db <- tryCatch(tools::Rd_db("svyrcs"), error = function(e) NULL)
  skip_if(is.null(db) || length(db) == 0)
  txt <- vapply(db, function(rd) paste(as.character(rd), collapse = " "), character(1L))

  expect_false(any(grepl("exactly `?degf\\(design\\)`? when there is nothing to impute", txt)))
  expect_false(any(grepl("exact when nothing is imputed", txt, fixed = TRUE)))

  ## and the corrected statement is actually present, so the claim cannot be fixed by deletion
  expect_true(any(grepl("nu + 1", txt, fixed = TRUE)))

  ## the k(m-1) <= 4 fallback is no longer sold as conservative. "conservative" on its own is not
  ## the tell -- the design-df default is legitimately described that way -- so this pins the
  ## sentence that disclaims it.
  expect_true(any(grepl("should not be read as conservative", txt, fixed = TRUE)))
})
