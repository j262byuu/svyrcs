# Plot a survey restricted cubic spline curve

Draws the estimated exposure-response curve with its confidence band, a
line at the null, and a marker at the reference value. Ratio measures
(hazard, odds and rate ratios) are drawn on a logarithmic y axis, which
is the scale on which the confidence band is symmetric; mean differences
are drawn on a linear axis.

## Usage

``` r
# S3 method for class 'svyrcs'
autoplot(
  object,
  xlab = NULL,
  ylab = NULL,
  title = NULL,
  rug = FALSE,
  band_alpha = 0.2,
  colour = "#2C6E9B",
  ...
)

# S3 method for class 'svyrcs'
plot(x, ...)
```

## Arguments

- xlab, ylab:

  Axis labels. The defaults name the exposure and the effect measure.

- title:

  Plot title. `NULL` (default) leaves it empty; `TRUE` builds one from
  the fit.

- rug:

  Draw a rug of the observed exposure values along the bottom axis. Off
  by default.

- band_alpha:

  Opacity of the confidence band.

- colour:

  Colour for the curve and the band.

- ...:

  Ignored.

- x, object:

  An object from
  [`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md).

## Value

A `ggplot` object.

## Details

The result is an ordinary `ggplot` object, so it can be modified in the
usual way, for example
`plot(fit) + ggplot2::labs(x = "Body mass index")`.

## Examples

``` r
design <- survey::svydesign(
  id = ~psu, strata = ~strata, weights = ~weight,
  nest = TRUE, data = nhanes_bmi
)
fit <- svyrcs(tchol ~ rcs(bmi, 4) + age + sex, design = design)
plot(fit)

plot(fit, title = TRUE) + ggplot2::labs(x = "Body mass index (kg/m2)")

#> NULL
```
