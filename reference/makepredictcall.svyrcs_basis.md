# Reuse fitted knots when predicting on new data

[`stats::makepredictcall()`](https://rdrr.io/r/stats/makepredictcall.html)
method for the basis returned by
[`rcs()`](https://j262byuu.github.io/svyrcs/reference/rcs.md). It
rewrites the
[`rcs()`](https://j262byuu.github.io/svyrcs/reference/rcs.md) call
stored in a model's terms so that it carries the knot locations chosen
at fit time. Without it, `predict(fit, newdata = ...)` would re-derive
knots from the quantiles of `newdata`, silently evaluating a different
spline.

## Usage

``` r
# S3 method for class 'svyrcs_basis'
makepredictcall(var, call)
```

## Arguments

- var:

  The fitted basis, as produced by
  [`rcs()`](https://j262byuu.github.io/svyrcs/reference/rcs.md).

- call:

  The call to
  [`rcs()`](https://j262byuu.github.io/svyrcs/reference/rcs.md) recorded
  in the model's terms.

## Value

The call with the fitted knots substituted for whatever `knots` was
originally given.

## Examples

``` r
fit <- lm(bmi ~ rcs(age, 4), data = nhanes_bmi)
# the terms now carry explicit knots rather than "4"
attr(terms(fit), "predvars")
#> list(bmi, rcs(age, knots = c(19, 38, 57, 80)))
```
