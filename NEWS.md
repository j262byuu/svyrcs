# svyrcs 0.1.0

First release.

* `svyrcs()` fits a restricted cubic spline model to a complex survey design and returns the
  exposure-response curve, tests and reference value in a single object. Survey-weighted Cox
  (`svycoxph`) and generalised linear (`svyglm`) models are supported, giving hazard ratios, odds
  ratios, rate ratios or mean differences.
* `rcs()` builds a restricted cubic spline basis for use inside a model formula. Knots are stored on
  the basis and reused by `predict()` through a `makepredictcall()` method, so predictions can never
  be made on silently re-derived knots. The basis is numerically identical to `rms::rcs()`.
* `svyrcs_curve()` estimates a curve from any already-fitted survey model.
* Curves are estimated as contrasts against a reference value, so covariates cancel exactly and no
  covariate reference grid is needed. Confidence intervals use a *t* quantile on the survey degrees
  of freedom from `survey::degf()`.
* Reference values may be given directly or chosen as the survey-weighted median, mean or an
  arbitrary weighted quantile, or as the minimum- or maximum-risk point of the curve.
* Tests of overall association and of non-linearity are reported as design-based *F* tests (matching
  `survey::regTermTest()`) alongside the chi-square form.
* `print()`, `summary()`, `plot()`, `autoplot()` and `predict()` methods.
* `nhanes_bmi`, a public-domain NHANES extract with linked mortality follow-up, for examples and
  tests.
