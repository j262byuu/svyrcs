## Base graphics implementation of the curve plot.
##
## This exists so that ggplot2 can live in Suggests without plotting becoming conditional: ggplot2
## contributes 15 recursive dependencies that survey does not already need, six of which need
## compilation, and a hard dependency is also a package whose CRAN archival would archive this one.

## Group colours. hcl.colors() is in base R and reasonably colourblind-safe; no attempt is made to
## imitate ggplot2's hue wheel.
group_palette <- function(n) {
  if (n <= 1L) return("#2C6E9B")
  grDevices::hcl.colors(n, "Dark 3")
}

## One panel: empty plot, then ribbons, reference lines, curves, markers, rug, legend. The order
## matters -- ribbons are drawn first so they never cover a line.
base_panel <- function(dat, cols, object, ylim, xlab, ylab, main, legend, rug, logax, band_alpha,
                       rug_values) {
  graphics::plot(range(dat$x), ylim, type = "n", log = logax, bty = "l",
                 xlab = xlab, ylab = ylab, main = main, las = 1)

  slice <- function(level) if (nzchar(level)) dat[dat$group == level, ] else dat

  for (i in seq_along(cols)) {
    d <- slice(names(cols)[i])
    if (!nrow(d)) next
    graphics::polygon(c(d$x, rev(d$x)), c(d$conf.low, rev(d$conf.high)),
                      col = grDevices::adjustcolor(cols[[i]], alpha.f = band_alpha), border = NA)
  }

  graphics::abline(h = object$null, lty = 2, col = "grey40")
  graphics::abline(v = object$ref$value, lty = 3, col = "grey40")

  for (i in seq_along(cols)) {
    d <- slice(names(cols)[i])
    if (!nrow(d)) next
    graphics::lines(d$x, d$estimate, col = cols[[i]], lwd = 2)
    graphics::points(object$ref$value, object$null, col = cols[[i]], pch = 19, cex = 1.1)
  }

  if (isTRUE(rug) && length(rug_values)) {
    inr <- rug_values >= min(dat$x) & rug_values <= max(dat$x)
    if (any(inr)) {
      graphics::rug(rug_values[inr], col = grDevices::adjustcolor("grey30", alpha.f = 0.15))
    }
  }

  if (isTRUE(legend) && length(cols) > 1L) {
    graphics::legend("topright", legend = names(cols), col = unlist(cols), lwd = 2, bty = "n",
                     title = object$groups$var, inset = 0.02, cex = 0.9)
  }
  invisible(NULL)
}

base_plot_svyrcs <- function(object, xlab = NULL, ylab = NULL, title = NULL, rug = FALSE,
                             band_alpha = 0.2, colour = "#2C6E9B", facet = FALSE) {
  cv <- as.data.frame(object$curve)
  ratio <- isTRUE(object$exponentiate)
  grouped <- !is.null(object$groups)

  xlab <- xlab %||% object$var
  ylab <- ylab %||% sprintf("%s (%g%% CI)", object$measure, 100 * object$level)
  logax <- if (ratio) "y" else ""

  levs <- if (grouped) object$groups$levels else ""
  cols <- if (grouped) as.list(group_palette(length(levs))) else list(colour)
  names(cols) <- levs

  ## One y range for every panel. Computing it per panel would give each subgroup its own axis,
  ## which makes a smaller effect look identical to a larger one -- the opposite of what panelling
  ## subgroups is for. ggplot2's facet_wrap() shares scales by default for the same reason.
  ylim <- range(c(cv$conf.low, cv$conf.high), finite = TRUE)
  if (ratio) ylim[1L] <- max(ylim[1L], .Machine$double.eps)

  rug_values <- if (isTRUE(rug)) {
    tryCatch(exposure_values(object$model, object$var, as_design_list(fit_design(object$model))),
             error = function(e) numeric(0))
  } else {
    numeric(0)
  }

  main <- if (isTRUE(title)) {
    sprintf("%s by %s", deparse1(object$formula[[2L]]), object$var)
  } else if (is.character(title)) {
    title
  } else {
    NULL
  }

  if (grouped && isTRUE(facet)) {
    old <- graphics::par(mfrow = c(1L, length(levs)), mar = c(4.2, 4.2, 2.4, 0.8))
    on.exit(graphics::par(old), add = TRUE)
    for (g in levs) {
      base_panel(cv[cv$group == g, ], cols[g], object, ylim, xlab, ylab, main = g,
                 legend = FALSE, rug = rug, logax = logax, band_alpha = band_alpha,
                 rug_values = rug_values)
    }
  } else {
    old <- graphics::par(mar = c(4.2, 4.2, if (is.null(main)) 1.2 else 3, 0.8))
    on.exit(graphics::par(old), add = TRUE)
    base_panel(cv, cols, object, ylim, xlab, ylab, main = main, legend = grouped,
               rug = rug, logax = logax, band_alpha = band_alpha, rug_values = rug_values)
  }

  invisible(object)
}
