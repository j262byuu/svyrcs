grp_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- svyrcs(tchol ~ rcs(bmi, 4) * sex + age, design = nhanes_design())
    }
    cached
  }
})

test_that("print reports the modifier, the interaction tests and a per-group table", {
  out <- capture.output(print(grp_fit()))
  expect_true(any(grepl("Modifier   sex, 2 levels, reference Male", out)))
  expect_true(any(grepl("Interaction", out)))
  expect_true(any(grepl("Shape interaction", out)))
  expect_true(any(grepl("Per group", out)))
  expect_true(any(grepl("^  Male ", out)))
  expect_true(any(grepl("^  Female ", out)))
})

test_that("summary reports the curve shape once per group", {
  s <- summary(grp_fit())
  expect_named(s$features, c("Male", "Female"))
  expect_equal(s$features$Male$min$estimate,
               min(grp_fit()$curve$estimate[grp_fit()$curve$group == "Male"]))

  out <- capture.output(print(s))
  expect_true(any(grepl("Curve shape: sex = Male", out)))
  expect_true(any(grepl("Curve shape: sex = Female", out)))
})

test_that("significant ranges are computed within a group, not across the stack", {
  fit <- grp_fit()
  s <- summary(fit)
  cv <- as.data.frame(fit$curve)
  for (g in fit$groups$levels) {
    sub <- cv[cv$group == g, ]
    sig <- s$features[[g]]$significant
    flagged <- rep(FALSE, nrow(sub))
    for (i in seq_len(nrow(sig))) {
      flagged <- flagged | (sub$x >= sig$from[i] & sub$x <= sig$to[i])
    }
    expect_equal(flagged, sub$conf.low > fit$null | sub$conf.high < fit$null, info = g)
  }
})

test_that("predict returns every group, or the one asked for", {
  fit <- grp_fit()
  p <- predict(fit, x = c(20, 30))
  expect_equal(nrow(p), 4L)
  expect_equal(levels(p$group), fit$groups$levels)

  one <- predict(fit, x = c(20, 30), group = "Female")
  expect_equal(nrow(one), 2L)
  expect_equal(as.character(unique(one$group)), "Female")
  expect_equal(one$estimate, p$estimate[p$group == "Female"])

  expect_error(predict(fit, x = 20, group = "Neither"), class = "svyrcs_error")
})

test_that("predict matches the stored curve within each group", {
  fit <- grp_fit()
  cv <- as.data.frame(fit$curve)
  for (g in fit$groups$levels) {
    sub <- cv[cv$group == g, ]
    at <- sub$x[c(1L, 75L, 200L)]
    p <- predict(fit, x = at, group = g)
    expect_equal(p$estimate, sub$estimate[c(1L, 75L, 200L)], info = g)
  }
})

test_that("autoplot maps colour and fill to the group and labels the legend", {
  skip_if_not_installed("ggplot2")
  p <- ggplot2::autoplot(grp_fit())
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$colour, "sex")
  expect_equal(p$labels$fill, "sex")
  expect_equal(nrow(p$data), 400L)
  expect_silent(invisible(ggplot2::ggplot_build(p)))
})

test_that("facet = TRUE panels the groups and drops the legend", {
  skip_if_not_installed("ggplot2")
  p <- ggplot2::autoplot(grp_fit(), facet = TRUE)
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(p$theme$legend.position, "none")
  expect_silent(invisible(ggplot2::ggplot_build(p)))

  expect_null(ggplot2::autoplot(grp_fit())$facet$params$facets$group)
})

test_that("the grouped subtitle carries the interaction p-value", {
  skip_if_not_installed("ggplot2")
  p <- ggplot2::autoplot(grp_fit(), title = TRUE)
  expect_match(p$labels$subtitle, "p\\(interaction by sex\\)")
  expect_lt(max(nchar(strsplit(p$labels$subtitle, "\n", fixed = TRUE)[[1L]])), 90)

  ## and an ungrouped fit does not mention it
  expect_false(grepl("interaction", ggplot2::autoplot(gaussian_fit(), title = TRUE)$labels$subtitle))
})

test_that("the ungrouped plot is unchanged", {
  skip_if_not_installed("ggplot2")
  p <- ggplot2::autoplot(gaussian_fit())
  expect_null(p$labels$colour)
  expect_null(p$labels$fill)
  expect_equal(nrow(p$data), 200L)
})
