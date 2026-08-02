# Restricted cubic spline analysis of complex survey data

Fits a restricted cubic spline model to a complex survey design and
returns the exposure-response curve, its confidence band, tests of
overall association and of non-linearity, and the reference value the
curve is anchored to.

## Usage

``` r
svyrcs(
  formula,
  design,
  family = NULL,
  ref = "median",
  ref_prob = 0.5,
  weighted_knots = TRUE,
  n = 200,
  range = NULL,
  level = 0.95,
  ...
)
```

## Arguments

- formula:

  A model formula containing exactly one
  [`rcs()`](https://j262byuu.github.io/svyrcs/reference/rcs.md) term for
  the exposure, plus any covariates. A
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) response
  selects a survey-weighted Cox model; anything else is fitted with
  [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html).

- design:

  A survey design object, as built by
  [`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html).
  Subset it with [`subset()`](https://rdrr.io/r/base/subset.html) before
  passing it in; the degrees of freedom are taken from the design you
  pass.

- family:

  Family for the
  [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html) path.
  Use [`quasibinomial()`](https://rdrr.io/r/stats/family.html) for a
  binary outcome (odds ratios),
  [`quasipoisson()`](https://rdrr.io/r/stats/family.html) for counts or
  rates (rate ratios), or the default
  [`gaussian()`](https://rdrr.io/r/stats/family.html) for a continuous
  outcome (mean differences). Ignored for Cox models.

- ref:

  Reference value the curve is anchored to: a number, or one of
  `"median"` (default), `"mean"`, `"quantile"`, `"min"` or `"max"`. See
  [references](https://j262byuu.github.io/svyrcs/reference/references.md).

- ref_prob:

  Probability used when `ref = "quantile"`.

- weighted_knots:

  If `TRUE` (default), a knot *count* is turned into knot locations
  using survey-weighted quantiles of the exposure, so knot placement
  reflects the population rather than the unweighted sample. Set to
  `FALSE` for unweighted quantiles. Explicit knot locations are always
  used as given.

- n:

  Number of points on the estimated curve.

- range:

  Length-2 numeric giving the exposure range to estimate over. Defaults
  to the 1st and 99th survey-weighted percentiles.

- level:

  Confidence level for the band.

- ...:

  Passed on to
  [`survey::svycoxph()`](https://rdrr.io/pkg/survey/man/svycoxph.html)
  or [`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html).

## Value

An object of class `svyrcs`, a list with components:

- curve:

  Data frame of `x`, `estimate`, `conf.low`, `conf.high`, `se`.

- ref:

  The reference value and how it was chosen.

- knots:

  Knot locations actually used.

- tests:

  Overall and non-linearity Wald tests, as *F* and chi-square.

- measure:

  `"HR"`, `"OR"`, `"RR"` or `"Difference"`.

- model:

  The fitted `svycoxph`/`svyglm` object, for your own diagnostics.

- degf, n, nevents, level, var, call:

  Analysis metadata.

## Why not just use `rms`

[`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html) gives the
same basis, but `rms` fitting functions do not know about sampling
weights, clustering or stratification, and `survey` has no spline
helper. `svyrcs` bridges the two: the spline is fitted inside a `survey`
model, so point estimates use the weights and the confidence band uses
Taylor linearisation on the design degrees of freedom.

## See also

[`svyrcs_curve()`](https://j262byuu.github.io/svyrcs/reference/svyrcs_curve.md)
for curves from a model you fitted yourself,
[`rcs()`](https://j262byuu.github.io/svyrcs/reference/rcs.md) for the
basis,
[references](https://j262byuu.github.io/svyrcs/reference/references.md)
for reference values.

## Examples

``` r
design <- survey::svydesign(
  id = ~psu, strata = ~strata, weights = ~weight,
  nest = TRUE, data = nhanes_bmi
)

# continuous outcome: mean difference in total cholesterol across the BMI range
fit <- svyrcs(tchol ~ rcs(bmi, 4) + age + sex, design = design)
fit
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

# survival outcome: all-cause mortality, anchored at the minimum-risk BMI
# \donttest{
library(survival)
fit_hr <- svyrcs(Surv(time, event) ~ rcs(bmi, 4) + age + sex,
                 design = design, ref = "min")
summary(fit_hr)
#> Restricted cubic spline on a complex survey design
#> 
#>   Model      survey-weighted Cox proportional hazards
#>   Outcome    Surv(time, event)
#>   Exposure   bmi, 4 knots at 19.78, 25.22, 29.71, 40.67 (weighted quantiles)
#>   Sample     10,617 observations, 1,893 events, 31 design df
#>   Reference  bmi = 27.13 (minimum-risk point), HR = 1
#> 
#>   Overall association  F = 48.07 on 3 and 31 df,  p = 9.1e-12 
#>   Non-linearity        F = 56.15 on 2 and 31 df,  p = 4.95e-11 
#> 
#>   Curve shape
#>     Lowest    bmi = 27.13,  HR = 1.00 (1.00, 1.00)
#>     Highest   bmi = 17.85,  HR = 2.85 (2.17, 3.75)
#>     95% band excludes HR = 1 over:
#>       bmi 17.85 to 25.40  (HR > 1)
#>       bmi 31.06 to 49.14  (HR > 1)
plot(fit_hr)


# binary outcome: odds of high total cholesterol
fit_or <- svyrcs(high_chol ~ rcs(bmi, 4) + age + sex, design = design,
                 family = quasibinomial())
# }
```
