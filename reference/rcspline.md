# Restricted cubic spline basis

Builds a restricted cubic spline (natural spline) basis for use inside a
model formula. The basis is Harrell's truncated power parameterisation
and matches [`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html) –
both the basis and, since 0.6.2, the knot placement this function
performs, including the small-sample rule that moves the outer knots to
the fifth smallest and fifth largest values below 100 observations.
Unlike `rms`, it stores its knots on the returned object and registers a
[`stats::makepredictcall()`](https://rdrr.io/r/stats/makepredictcall.html)
method, so a model fitted with `rcspline()` reuses exactly the same
knots when predicting on new data.

## Usage

``` r
rcspline(x, knots = 4)
```

## Arguments

- x:

  Numeric exposure vector.

- knots:

  Either the number of knots (3 to 7, default 4), placed at Harrell's
  recommended quantiles of `x`, or an explicit numeric vector of knot
  locations.
  [`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md)
  normally replaces a count with survey-weighted quantiles before
  fitting.

## Value

A numeric matrix with `length(x)` rows and `k - 1` columns, of class
`svyrcs_basis`, with attributes `knots` and `nk`.

## Details

Because a restricted cubic spline with `k` knots contributes `k - 1`
columns, the first of which is `x` itself, a test that the remaining
columns are jointly zero is a test of non-linearity. This is what
[`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md)
reports as the non-linearity *p*-value.

## Relationship to `rms`

This function used to be called `rcs()`, which masked
[`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html) whenever both
packages were attached. It was renamed in 0.7.0 because masking could
not be made safe: an `ols()`, `lrm()` or `cph()` model written with a
bare `rcs()` would pick up this basis, fit correctly, and then lose the
`Nonlinear` row from [`anova()`](https://rdrr.io/r/stats/anova.html)
with no error and no warning – and that row is usually the reason for
fitting a spline at all.

With the rename there is no collision. `rcs()` in an `rms` formula is
unambiguously
[`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html).

The basis itself is Harrell's, and agrees with
[`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html) to about
1e-14 **given the same knots**. Note the converse still holds: passing
`rcspline()` into an `rms` model fits, but the result is not an `rms`
design term, so [`anova()`](https://rdrr.io/r/stats/anova.html) reports
no `Nonlinear` row and `Predict()` fails. Use
[`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html) for `rms`
models and `rcspline()` for this package.

## See also

[`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md),
[`rcs_knots()`](https://j262byuu.github.io/svyrcs/reference/rcs_knots.md),
and [`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html) for `rms`
models – the two are not interchangeable in an `rms` formula, see the
section above.

## Examples

``` r
b <- svyrcs::rcspline(nhanes_bmi$bmi, 4)
dim(b)
#> [1] 10617     3
attr(b, "knots")
#> [1] 19.880 25.466 29.970 40.934
```
