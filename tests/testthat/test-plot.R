test_that("autoplot returns a ggplot whose data is the estimated curve", {
  fit <- gaussian_fit()
  p <- autoplot(fit)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data$x, fit$curve$x)
  expect_equal(p$data$estimate, fit$curve$estimate)
  expect_equal(p$data$conf.low, fit$curve$conf.low)
})

test_that("ratio measures get a log axis and differences do not", {
  skip_if_not_installed("survival")
  p_ratio <- autoplot(cox_fit())
  p_diff <- autoplot(gaussian_fit())
  has_log <- function(p) {
    any(vapply(p$scales$scales, function(s) {
      tr <- tryCatch(s$trans$name, error = function(e) NULL)
      isTRUE(grepl("log", tr %||% ""))
    }, logical(1L)))
  }
  expect_true(has_log(p_ratio))
  expect_false(has_log(p_diff))
})

test_that("axis labels default to the exposure and the measure, and can be overridden", {
  fit <- gaussian_fit()
  p <- autoplot(fit)
  expect_equal(p$labels$x, "bmi")
  expect_match(p$labels$y, "Difference")
  expect_match(p$labels$y, "95%")

  p2 <- autoplot(fit, xlab = "Body mass index", ylab = "Change in cholesterol")
  expect_equal(p2$labels$x, "Body mass index")
  expect_equal(p2$labels$y, "Change in cholesterol")
})

test_that("title = TRUE builds a wrapped subtitle, and a string sets the title", {
  fit <- gaussian_fit()
  expect_null(autoplot(fit)$labels$title)

  p <- autoplot(fit, title = TRUE)
  expect_match(p$labels$title, "bmi")
  expect_match(p$labels$subtitle, "design df")
  ## wrapped, so no single line is absurdly long
  expect_lt(max(nchar(strsplit(p$labels$subtitle, "\n", fixed = TRUE)[[1L]])), 90)

  expect_equal(autoplot(fit, title = "Custom")$labels$title, "Custom")
})

test_that("the plot builds without error, including the rug", {
  fit <- gaussian_fit()
  expect_silent(invisible(ggplot2::ggplot_build(autoplot(fit))))
  expect_silent(invisible(ggplot2::ggplot_build(autoplot(fit, rug = TRUE, title = TRUE))))
})

test_that("plot() draws and returns its argument invisibly", {
  fit <- gaussian_fit()
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  expect_invisible(plot(fit))
})
