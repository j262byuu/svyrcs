# Exposure-response curve from a fitted survey model

Estimates the exposure-response curve implied by an
[`rcspline()`](https://j262byuu.github.io/svyrcs/reference/rcspline.md)
term in an already-fitted survey model, as a contrast against a
reference exposure value.
[`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md)
calls this internally; use it directly when you have fitted the model
yourself, for instance because you need a model specification
[`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md) does
not build for you.

## Usage

``` r
svyrcs_curve(
  fit,
  var = NULL,
  ref = "median",
  ref_prob = 0.5,
  at = NULL,
  n = 200,
  range = NULL,
  level = 0.95,
  design = NULL,
  degf = NULL,
  group = NULL,
  df = NULL
)
```

## Arguments

- fit:

  A fitted model containing an
  [`rcspline()`](https://j262byuu.github.io/svyrcs/reference/rcspline.md)
  term, typically from
  [`survey::svycoxph()`](https://rdrr.io/pkg/survey/man/svycoxph.html)
  or [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html).
  Any model with [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`terms()`](https://rdrr.io/r/stats/terms.html) methods works.

- var:

  Name of the exposure. Only needed when the model contains more than
  one
  [`rcspline()`](https://j262byuu.github.io/svyrcs/reference/rcspline.md)
  term.

- ref:

  Reference value: a number, or one of `"median"`, `"mean"`,
  `"quantile"`, `"min"` (minimum-risk point) or `"max"` (maximum-risk
  point). See
  [`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md).

- ref_prob:

  Probability used when `ref = "quantile"`.

- at:

  Optional numeric vector of exposure values to evaluate. When given,
  `n` and `range` are ignored and the result has one row per element of
  `at`.

- n:

  Number of points on the curve.

- range:

  Length-2 numeric giving the exposure range to plot over. Defaults to
  the 1st and 99th survey-weighted percentiles of the exposure.

- level:

  Confidence level.

- design:

  The survey design. Defaults to the design stored on `fit`, which is
  what [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html)
  and
  [`survey::svycoxph()`](https://rdrr.io/pkg/survey/man/svycoxph.html)
  keep. Needed for weighted reference values and ranges.

- degf:

  Design degrees of freedom for the *t* quantile. Defaults to
  `survey::degf(design)`.

- group:

  When the spline is interacted with an effect modifier, the level or
  levels to estimate. `NULL` (default) returns every level, stacked,
  with a `group` column.

- df:

  Denominator degrees of freedom to use instead of the design's. `NULL`
  (default) keeps `survey::degf(design)`; a number substitutes it; `Inf`
  gives normal-quantile intervals. See the note in
  [`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md) on
  why this is a policy rather than a fixed answer.

## Value

A data frame of class `svyrcs_curve` with columns `x`, `estimate`,
`conf.low`, `conf.high` and `se` (on the linear-predictor scale), plus
`group` when the model has an effect modifier. Carries attributes `ref`,
`ref_method`, `measure`, `null`, `knots`, `var`, `degf`, `level` and
`modifier`.

## Details

The curve is computed from the spline coefficients alone. Every
covariate takes the same value at `x` and at the reference, so covariate
contributions cancel exactly and no covariate reference grid is needed.
Confidence intervals use a *t* quantile on the survey degrees of
freedom.

## See also

[`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md)

## Examples

``` r
design <- survey::svydesign(
  id = ~psu, strata = ~strata, weights = ~weight,
  nest = TRUE, data = nhanes_bmi
)
fit <- survey::svyglm(tchol ~ svyrcs::rcspline(bmi, 4) + age + sex, design = design)
curve <- svyrcs_curve(fit, "bmi", n = 20)
head(curve)
#>          x   estimate   conf.low   conf.high        se
#> 1 17.93000 -20.705258 -24.986868 -16.4236484 2.0993291
#> 2 19.56579 -16.030694 -19.215824 -12.8455636 1.5617109
#> 3 21.20158 -11.391725 -13.557005  -9.2264451 1.0616648
#> 4 22.83737  -7.080432  -8.460639  -5.7002245 0.6767336
#> 5 24.47316  -3.501327  -4.419794  -2.5828609 0.4503360
#> 6 26.10895  -1.051374  -1.549117  -0.5536305 0.2440500
```
