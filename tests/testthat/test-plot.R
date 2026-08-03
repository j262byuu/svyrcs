test_that("autoplot returns a ggplot whose data is the estimated curve", {
  skip_if_not_installed("ggplot2")
  fit <- gaussian_fit()
  p <- ggplot2::autoplot(fit)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data$x, fit$curve$x)
  expect_equal(p$data$estimate, fit$curve$estimate)
  expect_equal(p$data$conf.low, fit$curve$conf.low)
})

test_that("ratio measures get a log axis and differences do not", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")
  p_ratio <- ggplot2::autoplot(cox_fit())
  p_diff <- ggplot2::autoplot(gaussian_fit())
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
  skip_if_not_installed("ggplot2")
  fit <- gaussian_fit()
  p <- ggplot2::autoplot(fit)
  expect_equal(p$labels$x, "bmi")
  expect_match(p$labels$y, "Difference")
  expect_match(p$labels$y, "95%")

  p2 <- ggplot2::autoplot(fit, xlab = "Body mass index", ylab = "Change in cholesterol")
  expect_equal(p2$labels$x, "Body mass index")
  expect_equal(p2$labels$y, "Change in cholesterol")
})

test_that("title = TRUE builds a wrapped subtitle, and a string sets the title", {
  skip_if_not_installed("ggplot2")
  fit <- gaussian_fit()
  expect_null(ggplot2::autoplot(fit)$labels$title)

  p <- ggplot2::autoplot(fit, title = TRUE)
  expect_match(p$labels$title, "bmi")
  expect_match(p$labels$subtitle, "design df")
  ## wrapped, so no single line is absurdly long
  expect_lt(max(nchar(strsplit(p$labels$subtitle, "\n", fixed = TRUE)[[1L]])), 90)

  expect_equal(ggplot2::autoplot(fit, title = "Custom")$labels$title, "Custom")
})

test_that("the plot builds without error, including the rug", {
  skip_if_not_installed("ggplot2")
  fit <- gaussian_fit()
  expect_silent(invisible(ggplot2::ggplot_build(ggplot2::autoplot(fit))))
  expect_silent(invisible(ggplot2::ggplot_build(ggplot2::autoplot(fit, rug = TRUE, title = TRUE))))
})

test_that("plot() draws, and returns the ggplot invisibly so that `+` works", {
  skip_if_not_installed("ggplot2")
  fit <- gaussian_fit()
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)

  expect_invisible(plot(fit))
  p <- plot(fit)
  expect_s3_class(p, "ggplot")

  ## the documented idiom must actually modify the plot rather than evaluating to NULL, which is
  ## what happened while plot() returned the fit: ggplot2's `+.gg` swallows a non-ggplot silently
  modified <- plot(fit) + ggplot2::labs(x = "Body mass index")
  expect_s3_class(modified, "ggplot")
  expect_equal(modified$labels$x, "Body mass index")
})
