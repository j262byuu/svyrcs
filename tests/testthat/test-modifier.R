## `int_model()` lives in helper-fixture.R; alias it for brevity here.
int_fit <- function(rhs, data_extra = NULL) int_model(rhs, data_extra = data_extra)

test_that("no interaction means no modifier", {
  fit <- bare_glm()
  expect_null(find_modifier(fit, find_rcs_term(fit)))
})

test_that("a two-level factor modifier is found with its levels and columns", {
  fit <- int_fit("rcspline(bmi, knots = KN) * sex + age")
  term <- find_rcs_term(fit)
  m <- find_modifier(fit, term)

  expect_equal(m$var, "sex")
  expect_equal(m$levels, c("Male", "Female"))
  expect_equal(m$ref_level, "Male")
  expect_named(m$columns, "Female")
  expect_length(m$columns$Female, term$nk - 1L)
  expect_true(all(m$columns$Female %in% names(coef(fit))))
})

test_that("both formula orders find the same interaction columns", {
  a <- int_fit("rcspline(bmi, knots = KN) * sex + age")
  b <- int_fit("sex * rcspline(bmi, knots = KN) + age")

  ma <- find_modifier(a, find_rcs_term(a))
  mb <- find_modifier(b, find_rcs_term(b))

  expect_equal(ma$var, mb$var)
  expect_equal(ma$levels, mb$levels)
  ## the names themselves differ by orientation, which is exactly what has to be handled
  expect_true(all(grepl("^rcs", ma$columns$Female)))
  expect_true(all(grepl("^sexFemale", mb$columns$Female)))
  expect_true(all(mb$columns$Female %in% names(coef(b))))
})

test_that("a four-level factor gives one column set per non-reference level", {
  fit <- int_fit("rcspline(bmi, knots = KN) * race + age")
  term <- find_rcs_term(fit)
  m <- find_modifier(fit, term)

  expect_length(m$levels, 4L)
  expect_length(m$columns, 3L)
  expect_equal(names(m$columns), m$levels[-1L])
  for (cols in m$columns) {
    expect_length(cols, term$nk - 1L)
    expect_true(all(cols %in% names(coef(fit))))
  }
})

test_that("a logical modifier works", {
  fit <- int_fit("rcspline(bmi, knots = KN) * old + age",
                 data_extra = list(old = nhanes_bmi$age >= 65))
  m <- find_modifier(fit, find_rcs_term(fit))
  expect_equal(m$levels, c("FALSE", "TRUE"))
  expect_named(m$columns, "TRUE")
  expect_true(all(m$columns[["TRUE"]] %in% names(coef(fit))))
})

test_that("a numeric modifier is rejected with an explanation", {
  fit <- int_fit("rcspline(bmi, knots = KN) * age")
  expect_error(find_modifier(fit, find_rcs_term(fit)), class = "svyrcs_error")
  expect_error(find_modifier(fit, find_rcs_term(fit)), "Continuous effect modification")
})

test_that("a three-way interaction is rejected", {
  fit <- int_fit("rcspline(bmi, knots = KN) * sex * race")
  expect_error(find_modifier(fit, find_rcs_term(fit)), class = "svyrcs_error")
})

test_that("an interaction without the spline main effects is rejected", {
  fit <- int_fit("rcspline(bmi, knots = KN):sex + sex + age")
  expect_error(find_modifier(fit, find_rcs_term(fit)), class = "svyrcs_error")
  expect_error(find_modifier(fit, find_rcs_term(fit)), "main effect")
})

test_that("the selection matrix picks main effects, and adds interactions off the reference", {
  fit <- int_fit("rcspline(bmi, knots = KN) * sex + age")
  term <- find_rcs_term(fit)
  m <- find_modifier(fit, term)
  cn <- names(coef(fit))

  L_ref <- group_selection(cn, term, m, "Male")
  expect_equal(dim(L_ref), c(term$nk - 1L, length(cn)))
  expect_equal(rowSums(L_ref), setNames(rep(1, term$nk - 1L), rownames(L_ref)))

  L_oth <- group_selection(cn, term, m, "Female")
  expect_equal(rowSums(L_oth), setNames(rep(2, term$nk - 1L), rownames(L_oth)))
  expect_equal(which(L_oth[1L, ] == 1),
               setNames(match(c(m$spline_cols[1L], m$columns$Female[1L]), cn),
                        c(m$spline_cols[1L], m$columns$Female[1L])))

  ## reference-level selection equals the ungrouped selection
  expect_equal(L_ref, group_selection(cn, term, NULL, NULL))
})

test_that("the reference group's effective coefficients are the main effects", {
  fit <- int_fit("rcspline(bmi, knots = KN) * sex + age")
  term <- find_rcs_term(fit)
  m <- find_modifier(fit, term)
  b <- coef(fit)
  L <- group_selection(names(b), term, m, "Male")
  expect_equal(as.numeric(L %*% b), unname(b[m$spline_cols]))

  L2 <- group_selection(names(b), term, m, "Female")
  expect_equal(as.numeric(L2 %*% b),
               unname(b[m$spline_cols] + b[m$columns$Female]))
})
