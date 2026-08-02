# Changelog

## svyrcs 0.1.0

First release.

- [`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md)
  fits a restricted cubic spline model to a complex survey design and
  returns the exposure-response curve, tests and reference value in a
  single object. Survey-weighted Cox (`svycoxph`) and generalised linear
  (`svyglm`) models are supported, giving hazard ratios, odds ratios,
  rate ratios or mean differences.
- [`rcs()`](https://j262byuu.github.io/svyrcs/reference/rcs.md) builds a
  restricted cubic spline basis for use inside a model formula. Knots
  are stored on the basis and reused by
  [`predict()`](https://rdrr.io/r/stats/predict.html) through a
  [`makepredictcall()`](https://rdrr.io/r/stats/makepredictcall.html)
  method, so predictions can never be made on silently re-derived knots.
  The basis is numerically identical to
  [`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html).
- [`svyrcs_curve()`](https://j262byuu.github.io/svyrcs/reference/svyrcs_curve.md)
  estimates a curve from any already-fitted survey model.
- Curves are estimated as contrasts against a reference value, so
  covariates cancel exactly and no covariate reference grid is needed.
  Confidence intervals use a *t* quantile on the survey degrees of
  freedom from
  [`survey::degf()`](https://rdrr.io/pkg/survey/man/svychisq.html).
- Reference values may be given directly or chosen as the
  survey-weighted median, mean or an arbitrary weighted quantile, or as
  the minimum- or maximum-risk point of the curve.
- Tests of overall association and of non-linearity are reported as
  design-based *F* tests (matching
  [`survey::regTermTest()`](https://rdrr.io/pkg/survey/man/regTermTest.html))
  alongside the chi-square form.
- [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  and [`predict()`](https://rdrr.io/r/stats/predict.html) methods.
- `nhanes_bmi`, a public-domain NHANES extract with linked mortality
  follow-up, for examples and tests.
