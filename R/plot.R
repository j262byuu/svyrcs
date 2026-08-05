#' Plot a survey restricted cubic spline curve
#'
#' Draws the estimated exposure-response curve with its confidence band, a line at the null, and a
#' marker at the reference value. Ratio measures (hazard, odds and rate ratios) are drawn on a
#' logarithmic y axis, which is the scale on which the confidence band is symmetric; mean differences
#' are drawn on a linear axis.
#'
#' ggplot2 is used when it is installed, and an equivalent base graphics plot is drawn when it is
#' not, so plotting always works even though ggplot2 is only a suggested dependency. The two look
#' similar but are not pixel-identical; use `backend` to force one or the other.
#'
#' On the ggplot2 path the result is an ordinary `ggplot` object returned invisibly, so it can be
#' modified in the usual way: `plot(fit) + ggplot2::labs(x = "Body mass index")`.
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
#' @return
#' `autoplot()` always returns a `ggplot` object.
#'
#' `plot()` draws, and what it returns depends on the backend it used: the `ggplot` object
#' (invisibly) on the ggplot2 path, and the fit itself (invisibly) on the base graphics path. Code
#' that relies on the returned `ggplot`, such as `plot(fit) + ggplot2::labs(...)`, therefore needs
#' ggplot2 to be installed. Pass `backend = "ggplot2"` to make that requirement explicit and get a
#' clear error rather than a surprise.
#'
#' @examples
#' design <- survey::svydesign(
#'   id = ~psu, strata = ~strata, weights = ~weight,
#'   nest = TRUE, data = nhanes_bmi
#' )
#' # n = 50 rather than the default 200: the example draws the same curve and keeps well inside
#' # CRAN's five-second budget for a single example.
#' fit <- svyrcs(tchol ~ rcs(bmi, 4) + age + sex, design = design, n = 50)
#'
#' # works with or without ggplot2 installed
#' plot(fit)
#'
#' \donttest{
#' # the base fallback, and modifying the ggplot2 version. Both are cheap on their own; together
#' # with the fit above they push this one example past CRAN's five-second CPU budget, which is a
#' # limit per example rather than per file.
#' plot(fit, backend = "base", title = TRUE)
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot(fit, title = TRUE) + ggplot2::labs(x = "Body mass index (kg/m2)")
#' }
#' }
#'
#' @name plot.svyrcs
#'
#' @section Requires ggplot2:
#' `autoplot()` needs ggplot2, which is a suggested rather than a required dependency. Because of
#' that this package cannot re-export the generic, so `autoplot(fit)` needs `library(ggplot2)` first.
#' `plot(fit)` works either way, falling back to base graphics.
#'
#' @exportS3Method ggplot2::autoplot
autoplot.svyrcs <- function(object, xlab = NULL, ylab = NULL, title = NULL, rug = FALSE,
                            band_alpha = 0.2, colour = "#2C6E9B", facet = FALSE, ...) {
  if (!has_ggplot2()) {
    stop_svyrcs("autoplot() needs the 'ggplot2' package; install it with ",
                "install.packages(\"ggplot2\"), or use plot(), which falls back to base graphics")
  }
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

#' @param backend Which graphics system to draw with. `"auto"` (default) uses ggplot2 when it is
#'   installed and base graphics otherwise; `"ggplot2"` and `"base"` force one or the other. Forcing
#'   `"ggplot2"` without the package installed is an error.
#'
#' @rdname plot.svyrcs
#' @export
plot.svyrcs <- function(x, ..., backend = c("auto", "ggplot2", "base")) {
  backend <- match.arg(backend)
  have_ggplot2 <- has_ggplot2()

  if (backend == "ggplot2" && !have_ggplot2) {
    stop_svyrcs("`backend = \"ggplot2\"` needs the 'ggplot2' package; install it with ",
                "install.packages(\"ggplot2\"), or use backend = \"base\"")
  }

  if (backend == "base" || (backend == "auto" && !have_ggplot2)) {
    return(base_plot_svyrcs(x, ...))
  }

  ## Return the ggplot, not the fit. Returning the fit is what `plot(fit) + ggplot2::labs(...)`
  ## needs, and returning the fit instead made that idiom evaluate to NULL: ggplot2's `+.gg` accepts
  ## a non-ggplot left operand and silently yields NULL rather than erroring, so neither the user nor
  ## R CMD check ever saw a problem. Invisibly, so `plot(fit)` at the console still just draws.
  p <- autoplot.svyrcs(x, ...)
  print(p)
  invisible(p)
}
