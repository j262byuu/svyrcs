# Reuse fitted knots when predicting on new data

[`stats::makepredictcall()`](https://rdrr.io/r/stats/makepredictcall.html)
method for the basis returned by
[`rcspline()`](https://j262byuu.github.io/svyrcs/reference/rcspline.md).
It rewrites the
[`rcspline()`](https://j262byuu.github.io/svyrcs/reference/rcspline.md)
call stored in a model's terms so that it carries the knot locations
chosen at fit time. Without it, `predict(fit, newdata = ...)` would
re-derive knots from the quantiles of `newdata`, silently evaluating a
different spline.

## Usage

``` r
# S3 method for class 'svyrcs_basis'
makepredictcall(var, call)
```

## Arguments

- var:

  The fitted basis, as produced by
  [`rcspline()`](https://j262byuu.github.io/svyrcs/reference/rcspline.md).

- call:

  The call to
  [`rcspline()`](https://j262byuu.github.io/svyrcs/reference/rcspline.md)
  recorded in the model's terms.

## Value

The call with the fitted knots substituted for whatever `knots` was
originally given.

## Examples

``` r
fit <- lm(bmi ~ svyrcs::rcspline(age, 4), data = nhanes_bmi)
# the terms now carry explicit knots rather than "4"
attr(terms(fit), "predvars")
#> list(bmi, svyrcs::rcspline(age, knots = c(19, 38, 57, 80)))
```
