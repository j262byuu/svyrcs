#' Plot a survey restricted cubic spline curve
#'
#' Draws the estimated exposure-response curve with its confidence band, a line at the null, and a
#' marker at the reference value. Ratio measures (hazard, odds and rate ratios) are drawn on a
#' logarithmic y axis, which is the scale on which the confidence band is symmetric; mean differences
#' are drawn on a linear axis.
#'
#' The result is an ordinary `ggplot` object, so it can be modified in the usual way, for example
#' `plot(fit) + ggplot2::labs(x = "Body mass index")`.
#'
#' @param x,object An object from [svyrcs()].
#' @param xlab,ylab Axis labels. The defaults name the exposure and the effect measure.
#' @param title Plot title. `NULL` (default) leaves it empty; `TRUE` builds one from the fit.
#' @param rug Draw a rug of the observed exposure values along the bottom axis. Off by default.
#' @param band_alpha Opacity of the confidence band.
#' @param colour Colour for the curve and the band. Ignored when the fit has an effect modifier, in
#'   which case groups are coloured by the default discrete scale.
#' @param facet For a fit with an effect modifier, `TRUE` puts each group in its own panel instead of
#'   overlaying them.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' design <- survey::svydesign(
#'   id = ~psu, strata = ~strata, weights = ~weight,
#'   nest = TRUE, data = nhanes_bmi
#' )
#' fit <- svyrcs(tchol ~ rcs(bmi, 4) + age + sex, design = design)
#' plot(fit)
#' plot(fit, title = TRUE) + ggplot2::labs(x = "Body mass index (kg/m2)")
#'
#' @name plot.svyrcs
#' @export
autoplot.svyrcs <- function(object, xlab = NULL, ylab = NULL, title = NULL, rug = FALSE,
                            band_alpha = 0.2, colour = "#2C6E9B", facet = FALSE, ...) {
  cv <- as.data.frame(object$curve)
  ratio <- isTRUE(object$exponentiate)
  grouped <- !is.null(object$groups)

  xlab <- xlab %||% object$var
  ylab <- ylab %||% sprintf("%s (%g%% CI)", object$measure, 100 * object$level)

  ## The reference marker is one point per panel when facetting, and a single shared point
  ## otherwise, so it has to carry the grouping variable in the facetted case.
  ref_pt <- data.frame(x = object$ref$value, estimate = object$null)
  if (grouped) {
    ref_pt <- do.call(rbind, lapply(object$groups$levels, function(g) {
      cbind(ref_pt, group = factor(g, levels = object$groups$levels))
    }))
  }

  p <- ggplot2::ggplot(cv, ggplot2::aes(x = .data$x, y = .data$estimate))

  if (grouped) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high, fill = .data$group),
        alpha = band_alpha, colour = NA
      ) +
      ggplot2::geom_hline(yintercept = object$null, linetype = "dashed",
                          colour = "grey40", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = object$ref$value, linetype = "dotted",
                          colour = "grey40", linewidth = 0.4) +
      ggplot2::geom_line(ggplot2::aes(colour = .data$group), linewidth = 0.9) +
      ggplot2::geom_point(data = ref_pt, ggplot2::aes(colour = .data$group), size = 2.2,
                          show.legend = FALSE) +
      ggplot2::labs(colour = object$groups$var, fill = object$groups$var)
  } else {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high),
                           fill = colour, alpha = band_alpha) +
      ggplot2::geom_hline(yintercept = object$null, linetype = "dashed",
                          colour = "grey40", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = object$ref$value, linetype = "dotted",
                          colour = "grey40", linewidth = 0.4) +
      ggplot2::geom_line(colour = colour, linewidth = 0.9) +
      ggplot2::geom_point(data = ref_pt, colour = colour, size = 2.2)
  }

  p <- p +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.4),
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "grey30"),
      strip.background = ggplot2::element_blank(),
      legend.position = if (grouped && !isTRUE(facet)) "right" else "none"
    )

  if (grouped && isTRUE(facet)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$group))
  }

  if (ratio) {
    p <- p + ggplot2::scale_y_continuous(transform = "log10")
  }

  if (isTRUE(rug)) {
    xv <- tryCatch(exposure_values(object$model, object$var, fit_design(object$model)),
                   error = function(e) NULL)
    if (!is.null(xv)) {
      inr <- xv >= min(cv$x) & xv <= max(cv$x)
      p <- p + ggplot2::geom_rug(data = data.frame(x = xv[inr], estimate = object$null),
                                 sides = "b", alpha = 0.08, colour = "grey30")
    }
  }

  if (isTRUE(title)) {
    ## Wrapped, because an unwrapped subtitle runs off the edge of the panel and is silently
    ## clipped -- which looks like a rendering bug rather than a long string.
    subtitle <- sprintf(
      "%s, %d knots, %d design df; reference %s = %s (%s); p(overall) %s, p(non-linearity) %s%s",
      model_description(object), object$nk, round(object$degf), object$var,
      fmt_num(object$ref$value, 4), object$ref$method,
      fmt_p(object$tests$overall$p_F), fmt_p(object$tests$nonlinear$p_F),
      ## For a grouped fit the interaction p-value is the number a reader looks for first.
      if (grouped) {
        sprintf("; p(interaction by %s) %s", object$groups$var,
                fmt_p(object$tests$interaction$p_F))
      } else {
        ""
      }
    )
    p <- p + ggplot2::labs(
      title = sprintf("%s by %s", deparse1(object$formula[[2L]]), object$var),
      subtitle = paste(strwrap(subtitle, width = 72), collapse = "\n")
    )
  } else if (is.character(title)) {
    p <- p + ggplot2::labs(title = title)
  }

  p
}

#' @rdname plot.svyrcs
#' @export
plot.svyrcs <- function(x, ...) {
  print(autoplot(x, ...))
  invisible(x)
}

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
