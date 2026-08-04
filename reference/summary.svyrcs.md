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
two tests and a `features` component holding the grid minimum, the grid
maximum, and `significant`: the stretches over which the **pointwise**
band excludes the null.

## These ranges are pointwise

The band is a pointwise confidence band: each grid point has its own
interval at the stated level. A stretch where those intervals exclude
the null is therefore **not** a simultaneous confidence region and
carries no family-wise error control. With 200 grid points, some will
cross a 95% pointwise band by chance even when the true curve is flat.
Read the ranges as a description of where the curve sits relative to the
null, not as a multiplicity-corrected finding.

The reported lowest and highest points are likewise the extremes of the
evaluated grid, so they can move a little if `n` changes.

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
#>   Exposure   bmi, 4 knots at 19.85, 25.27, 29.71, 40.65 (weighted quantiles)
#>   Sample     9,968 observations, 31 design df
#>   Reference  bmi = 27.37 (weighted median), Difference = 0
#>   Dropped    649 rows with missing values
#> 
#>   Overall association  F = 39.86 on 3 and 31 df,  p = 9.36e-11 
#>   Non-linearity        F = 55.21 on 2 and 31 df,  p = 6.07e-11 
#> 
#>   Curve shape
#>     Grid lowest   bmi = 17.93,  Difference = -20.7 (-25.0, -16.4)
#>     Grid highest  bmi = 30.11,  Difference = 0.670 (-0.587, 1.93)
#>     Pointwise 95% band excludes Difference = 0 over:
#>       bmi 17.93 to 27.30  (Difference < 0)
#>       bmi 27.46 to 28.08  (Difference > 0)
#>       bmi 38.08 to 49.01  (Difference < 0)
```
