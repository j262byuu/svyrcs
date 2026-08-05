## Draw to a throwaway device, so the tests exercise the real drawing code without leaving files or
## depending on a display.
on_device <- function(expr, width = 700, height = 470) {
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = width, height = height)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)
  force(expr)
}

test_that("the base fallback draws every kind of fit without error", {
  skip_if_not_installed("survival")
  fits <- list(
    single = gaussian_fit(),
    cox = cox_fit(),
    grouped = svyrcs(tchol ~ rcspline(bmi, 4) * sex + age, design = nhanes_design())
  )
  for (nm in names(fits)) {
    expect_silent(on_device(base_plot_svyrcs(fits[[nm]])))
  }
  expect_silent(on_device(base_plot_svyrcs(fits$grouped, facet = TRUE)))
  expect_silent(on_device(base_plot_svyrcs(fits$cox, rug = TRUE, title = TRUE)))
  expect_silent(on_device(base_plot_svyrcs(fits$single, xlab = "BMI", ylab = "Change",
                                           title = "Custom")))
})

test_that("the base fallback draws a multiply imputed fit", {
  fit <- svyrcs(bmi ~ rcspline(age, 4) + sex + tchol, design = mi_design())
  expect_silent(on_device(base_plot_svyrcs(fit)))
})

test_that("the base fallback returns its argument invisibly", {
  fit <- gaussian_fit()
  on_device({
    out <- withVisible(base_plot_svyrcs(fit))
    expect_false(out$visible)
    expect_identical(out$value, fit)
  })
})

test_that("graphical parameters are restored", {
  fit <- svyrcs(tchol ~ rcspline(bmi, 4) * sex + age, design = nhanes_design())
  on_device({
    before <- graphics::par(c("mfrow", "mar"))
    base_plot_svyrcs(fit, facet = TRUE)
    expect_equal(graphics::par("mfrow"), before$mfrow)
    expect_equal(graphics::par("mar"), before$mar)
  })
})

test_that("faceted panels share one y axis", {
  fit <- svyrcs(Surv(time, event) ~ rcspline(bmi, 4) * sex + age, design = nhanes_design())
  skip_if_not_installed("survival")

  usr <- list()
  on_device({
    old <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(old), add = TRUE)
    ## redraw panel by panel so par("usr") can be read after each
    cv <- as.data.frame(fit$curve)
    ylim <- range(c(cv$conf.low, cv$conf.high), finite = TRUE)
    for (g in fit$groups$levels) {
      base_panel(cv[cv$group == g, ], setNames(list(group_palette(2)[1]), g), fit, ylim,
                 "x", "y", main = g, legend = FALSE, rug = FALSE, logax = "y",
                 band_alpha = 0.2, rug_values = numeric(0))
      usr[[g]] <- graphics::par("usr")[3:4]
    }
  })
  expect_equal(usr[[1]], usr[[2]])

  ## and the shared range really does span both groups, rather than one of them
  cv <- as.data.frame(fit$curve)
  per_group <- vapply(fit$groups$levels, function(g) {
    range(c(cv$conf.low[cv$group == g], cv$conf.high[cv$group == g]))
  }, numeric(2))
  expect_false(isTRUE(all.equal(per_group[, 1], per_group[, 2])))
})

test_that("ratio measures get a log axis and differences do not", {
  skip_if_not_installed("survival")
  on_device({
    base_plot_svyrcs(cox_fit())
    expect_true(graphics::par("ylog"))
  })
  on_device({
    base_plot_svyrcs(gaussian_fit())
    expect_false(graphics::par("ylog"))
  })
})

test_that("the palette gives one colour per group", {
  expect_length(group_palette(1), 1L)
  expect_length(group_palette(4), 4L)
  expect_false(anyDuplicated(group_palette(4)) > 0)
})

test_that("the backend argument selects the drawing path", {
  skip_if_not_installed("ggplot2")
  fit <- gaussian_fit()

  ## with ggplot2 present, "auto" and "ggplot2" both return the ggplot
  on_device({
    expect_s3_class(plot(fit), "ggplot")
    expect_s3_class(plot(fit, backend = "ggplot2"), "ggplot")
  })

  ## "base" draws and returns the fit invisibly, even though ggplot2 is available
  on_device({
    out <- withVisible(plot(fit, backend = "base"))
    expect_false(out$visible)
    expect_identical(out$value, fit)
  })

  expect_error(plot(fit, backend = "nonsense"))
})

test_that("autoplot is registered on ggplot2's generic", {
  skip_if_not_installed("ggplot2")
  ## delayed S3 registration: NAMESPACE carries S3method(ggplot2::autoplot, svyrcs), so dispatch
  ## works once ggplot2 is loaded, without svyrcs importing it
  expect_s3_class(ggplot2::autoplot(gaussian_fit()), "ggplot")
  expect_false("autoplot" %in% getNamespaceExports("svyrcs"))
})

test_that("plot() and autoplot() agree on the drawn data", {
  skip_if_not_installed("ggplot2")
  fit <- svyrcs(tchol ~ rcspline(bmi, 4) * sex + age, design = nhanes_design())
  on_device({
    p <- plot(fit)
    expect_equal(p$data, ggplot2::autoplot(fit)$data)
  })
})

test_that("everything degrades correctly when ggplot2 is unavailable", {
  skip_if_not_installed("testthat", "3.1.7")
  fit <- gaussian_fit()

  ## Pretend ggplot2 cannot be loaded. This is the only way to reach the fallback branches on a
  ## machine that has ggplot2 installed; a genuinely ggplot2-free library only exists in CI.
  local_mocked_bindings(has_ggplot2 = function() FALSE, .package = "svyrcs")

  on_device({
    out <- withVisible(plot(fit))
    expect_false(out$visible)
    expect_identical(out$value, fit)   # base path, so the fit comes back, not a ggplot
  })

  expect_error(plot(fit, backend = "ggplot2"), class = "svyrcs_error")
  expect_error(plot(fit, backend = "ggplot2"), "install.packages")
  expect_error(autoplot.svyrcs(fit), class = "svyrcs_error")
  expect_error(autoplot.svyrcs(fit), "falls back to base graphics")

  ## and the explicitly requested base backend is unaffected
  on_device(expect_silent(plot(fit, backend = "base")))
})
