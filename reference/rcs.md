# Restricted cubic spline basis

Builds a restricted cubic spline (natural spline) basis for use inside a
model formula. The basis is Harrell's truncated power parameterisation
and is numerically identical to
[`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html), but it
stores its knots on the returned object and registers a
[`stats::makepredictcall()`](https://rdrr.io/r/stats/makepredictcall.html)
method, so a model fitted with `rcs()` reuses exactly the same knots
when predicting on new data.

## Usage

``` r
rcs(x, knots = 4)
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

## Masking `rms`

Attaching `svyrcs` after `rms` masks
[`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html). The two
bases agree to about 1e-14, so results are unchanged; use
[`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html) explicitly if
you need the `rms` version.

## See also

[`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md),
[`rcs_knots()`](https://j262byuu.github.io/svyrcs/reference/rcs_knots.md)

## Examples

``` r
b <- rcs(nhanes_bmi$bmi, 4)
dim(b)
#> [1] 10617     3
attr(b, "knots")
#> [1] 19.880 25.466 29.970 40.934
```
