# Effect estimates at specific exposure values

Gives the estimated effect, with confidence interval, at exposure values
you choose, rather than on the plotting grid. This is what you want when
a paper reports, say, the hazard ratio at a BMI of 30 against a BMI of
25.

## Usage

``` r
# S3 method for class 'svyrcs'
predict(object, x, ref = NULL, level = NULL, group = NULL, ...)
```

## Arguments

- object:

  An object from
  [`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md).

- x:

  Numeric vector of exposure values.

- ref:

  Reference value to contrast against. Defaults to the reference the fit
  was anchored to; supply a number to compare against a different value
  without refitting.

- level:

  Confidence level. Defaults to the level used by the fit.

- group:

  For a fit with an effect modifier, the level or levels wanted. `NULL`
  (default) returns every level.

- ...:

  Ignored.

## Value

A data frame with columns `x`, `estimate`, `conf.low`, `conf.high` and
`se`, plus `group` when the fit has an effect modifier.

## Examples

``` r
design <- survey::svydesign(
  id = ~psu, strata = ~strata, weights = ~weight,
  nest = TRUE, data = nhanes_bmi
)
fit <- svyrcs(tchol ~ svyrcs::rcspline(bmi, 4) + age + sex, design = design)

# against the fitted reference
predict(fit, x = c(20, 25, 30, 35))
#>    x    estimate    conf.low  conf.high        se
#> 1 20 -14.7250883 -17.6102495 -11.839927 1.4146321
#> 2 25  -2.5271192  -3.3727270  -1.681511 0.4146125
#> 3 30   0.6693980  -0.5452256   1.884022 0.5955458
#> 4 35  -0.6040603  -2.7862079   1.578087 1.0699354

# against a BMI of 25 instead
predict(fit, x = c(20, 30, 35), ref = 25)
#>    x   estimate    conf.low conf.high        se
#> 1 20 -12.197969 -14.9313651 -9.464573 1.3402196
#> 2 30   3.196517   1.1953994  5.197635 0.9811741
#> 3 35   1.923059  -0.9642341  4.810352 1.4156774
```
