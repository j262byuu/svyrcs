# svyrcs 0.2.0

* `svyrcs()` accepts an interaction between the spline and a grouping variable,
  `rcs(bmi, 4) * sex`, and returns one curve per group from a single model. Writing the terms the
  other way round, `sex * rcs(bmi, 4)`, works too.
* Two effect-modification tests are reported: **Interaction**, whether the association differs
  between groups at all, and **Shape interaction**, whether the curvature differs as opposed to the
  whole curve being shifted. The interaction test matches `survey::regTermTest()`. Overall and
  non-linearity tests are also reported within each group.
* Group curves are estimated through a selection matrix on the full coefficient vector, so the
  standard error for a non-reference group correctly includes the covariance between the main and
  interaction terms.
* All groups share one reference value so the curves are comparable. `ref = "min"` and `ref = "max"`
  are located on the reference level's curve, and the output says which level that was.
* `predict()` gains a `group` argument; `plot()` and `autoplot()` colour by group and gain
  `facet = TRUE`; `print()` and `summary()` report per-group results.
* A rank-deficient interaction — the spline collinear with the modifier inside a group, typically
  when the modifier is derived from the exposure — is now diagnosed as such rather than reported as
  a missing coefficient.
* Fits without an interaction are unchanged, and a regression test pins their numbers and printed
  output.

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
