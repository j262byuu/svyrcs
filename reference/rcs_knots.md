# Knot locations for a restricted cubic spline

Resolves the `knots` argument of
[`rcs()`](https://j262byuu.github.io/svyrcs/reference/rcs.md) into an
explicit, sorted vector of knot locations. Exported so that knot
placement can be inspected, or computed once and reused across models.

## Usage

``` r
rcs_knots(x, knots = 4, var = "x")
```

## Arguments

- x:

  Numeric vector of exposure values.

- knots:

  Either a single count between 3 and 7, in which case knots are placed
  at Harrell's recommended quantiles of `x`, or a numeric vector of at
  least three explicit knot locations.

- var:

  Optional variable name, used only to make error messages readable.

## Value

A sorted numeric vector of knot locations.

## Examples

``` r
rcs_knots(nhanes_bmi$bmi, 4)
#> [1] 19.880 25.466 29.970 40.934
rcs_knots(nhanes_bmi$bmi, c(20, 25, 30, 40))
#> [1] 20 25 30 40
```
