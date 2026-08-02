# Effect estimates at specific exposure values

Gives the estimated effect, with confidence interval, at exposure values
you choose, rather than on the plotting grid. This is what you want when
a paper reports, say, the hazard ratio at a BMI of 30 against a BMI of
25.

## Usage

``` r
# S3 method for class 'svyrcs'
predict(object, x, ref = NULL, level = NULL, ...)
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

- ...:

  Ignored.

## Value

A data frame with columns `x`, `estimate`, `conf.low`, `conf.high` and
`se`.

## Examples

``` r
design <- survey::svydesign(
  id = ~psu, strata = ~strata, weights = ~weight,
  nest = TRUE, data = nhanes_bmi
)
fit <- svyrcs(tchol ~ rcs(bmi, 4) + age + sex, design = design)

# against the fitted reference
predict(fit, x = c(20, 25, 30, 35))
#>    x    estimate    conf.low  conf.high        se
#> 1 20 -14.6941295 -17.5737038 -11.814555 1.4118928
#> 2 25  -2.5053629  -3.3460312  -1.664695 0.4121906
#> 3 30   0.6924554  -0.5263881   1.911299 0.5976148
#> 4 35  -0.5684705  -2.7505635   1.613622 1.0699086

# against a BMI of 25 instead
predict(fit, x = c(20, 30, 35), ref = 25)
#>    x   estimate    conf.low conf.high        se
#> 1 20 -12.188767 -14.9197540 -9.457779 1.3390387
#> 2 30   3.197818   1.1960692  5.199567 0.9814837
#> 3 35   1.936892  -0.9476258  4.821411 1.4143168
```
