# Summarise a survey restricted cubic spline fit

In addition to what [`print()`](https://rdrr.io/r/base/print.html)
shows, the summary reports the shape of the fitted curve: where the
estimated effect is smallest and largest, and over which stretches of
the exposure the confidence band excludes the null.

## Usage

``` r
# S3 method for class 'svyrcs'
summary(object, ...)
```

## Arguments

- object:

  An object from
  [`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md).

- ...:

  Ignored.

## Value

An object of class `summary.svyrcs`, a list with the fit metadata, the
two tests and a `features` component holding the minimum, the maximum
and the significant ranges.

## Examples

``` r
design <- survey::svydesign(
  id = ~psu, strata = ~strata, weights = ~weight,
  nest = TRUE, data = nhanes_bmi
)
fit <- svyrcs(tchol ~ rcs(bmi, 4) + age + sex, design = design)
summary(fit)
#> Restricted cubic spline on a complex survey design
#> 
#>   Model      survey-weighted GLM (gaussian, identity link)
#>   Outcome    tchol
#>   Exposure   bmi, 4 knots at 19.78, 25.22, 29.71, 40.67 (weighted quantiles)
#>   Sample     9,968 observations, 31 design df
#>   Reference  bmi = 27.35 (weighted median), Difference = 0
#>   Dropped    649 rows with missing values
#> 
#>   Overall association  F = 39.87 on 3 and 31 df,  p = 9.34e-11 
#>   Non-linearity        F = 55.21 on 2 and 31 df,  p = 6.07e-11 
#> 
#>   Curve shape
#>     Lowest    bmi = 17.85,  Difference = -20.9 (-25.2, -16.5)
#>     Highest   bmi = 30.11,  Difference = 0.694 (-0.568, 1.96)
#>     95% band excludes Difference = 0 over:
#>       bmi 17.85 to 27.28  (Difference < 0)
#>       bmi 27.44 to 28.07  (Difference > 0)
#>       bmi 38.13 to 49.14  (Difference < 0)
```
