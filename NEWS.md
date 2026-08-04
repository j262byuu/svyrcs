# svyrcs 0.6.0

## Multiple-imputation joint tests now use a named rule

**Joint-test *p*-values under imputation change.** They previously came from a construction of my
own: a Wald statistic on the pooled total covariance, with a denominator degrees of freedom obtained
by collapsing the block's relative increase in variance to a scalar and reusing the scalar
Barnard-Rubin formula. That is not D1, D2 or D3, and a reader could not tell what had been computed.

The overall, non-linearity, interaction and per-group tests now use **D1** (Li, Raghunathan and
Rubin 1991) with the small-sample denominator degrees of freedom of **Reiter (2007)**, which is the
correction designed for a finite complete-data degrees of freedom — the survey situation. The
implementation agrees with `mitml:::.D1` to machine precision, which the test suite checks.

Two differences from the old rule, both of which move results:

* the statistic uses \eqn{(1 + r)\bar{U}} rather than the pooled total covariance;
* the denominator degrees of freedom come from Reiter's expansion.

Where \eqn{k(m-1) \le 4} Reiter's expansion is undefined — it divides by \eqn{k(m-1) - 4}, and
`mitml` returns a degenerate value there. `svyrcs` warns and falls back to the complete-data degrees
of freedom; raising the number of imputations is the fix.

## Two related corrections

* **`df` now reaches the imputation arithmetic.** In 0.5.1 it affected only the single-design path,
  so `df = Inf` silently did nothing to a multiply imputed fit.
* **`barnard_rubin_df()` follows the published formula at zero between-imputation variance**, where
  it returns \eqn{(\nu+1)\nu/(\nu+3)} rather than \eqn{\nu}. Earlier versions short-circuited
  to \eqn{\nu} so that a fit with nothing to impute reproduced the complete-data fit exactly. That
  was a departure from the rule dressed up as an improvement. The consequence: with nothing to
  impute, estimates, standard errors and the *F* statistic are still identical to the complete-data
  fit, but intervals are about 1% wider. The scalar formula now matches `mice:::barnard.rubin` for
  every combination tested.

## What the simulation showed

A null simulation, kept with the development notes, found the imputation joint test **conservative**:
0.7% to 3.2% rejection at a nominal 5%, against 5.8% to 9.0% for the same test on complete data, with
pointwise coverage of 96.6% to 98.8% against 95%. It is worst with few primary sampling units. This
is now documented in `?svyrcs` and the vignette, because a non-significant imputation *p*-value is
weaker evidence than it looks.

This is not a reason to revert: the previous rule was equally conservative or more so on the same
data, and had no name to look up. The change buys a citable rule at no measured cost in size.

The first version of that simulation was itself wrong, and the way it failed is worth recording. It
reported every rule as conservative — including a naive chi-square, which cannot be. The exposure
had been made missing at random and imputed from the outcome, but under the null those are
independent, so the imputation model had no signal and the between-imputation variance swamped
everything. Running the same generator with the missingness switched off showed the harness was fine
and localised the fault.

The per-point curve degrees of freedom are unchanged in method — scalar Barnard-Rubin remains the
right rule for a scalar estimand.

# svyrcs 0.5.1

Labelling and documentation, from the same review as 0.5.0. **No estimate, interval or p-value
changes**, with one exception noted below; the compat suite passes unmodified.

The theme is labels that claimed more than the method behind them delivered.

* **Ranges in `summary()` are now called pointwise.** They come from a pointwise band, so a stretch
  where the intervals exclude the null is not a simultaneous region and has no family-wise error
  control — with 200 grid points, some will cross a 95% band by chance even when the curve is flat.
* **`ref = "min"` and `ref = "max"` state that the anchor is selected from the same data**, and that
  the uncertainty in its location is not propagated. At the selected point the estimate is the null
  with zero standard error by construction; that is an artefact of anchoring there, not evidence.
* **The turning points reported by `summary()` are labelled grid extrema**, because they are the
  extremes of the evaluated curve and can shift when `n` changes.
* **New `df` argument** on `svyrcs()` and `svyrcs_curve()`. The design degrees of freedom remain the
  default, but they are a conservative policy rather than the only correct answer, and the
  documentation now says so. `df = Inf` gives normal-quantile intervals. Under multiple imputation
  `df` supplies the complete-data degrees of freedom the Barnard-Rubin correction starts from.
* **Variance-method wording fixed.** DESCRIPTION, README, the manual and the vignette described
  everything as Taylor linearisation, but replicate-weight designs are accepted and their variance
  comes from the replicate weights.
* **A Poisson model without an offset now reports a mean ratio rather than a rate ratio** — the only
  labelling change that alters output. With person-time in an offset it reports a rate ratio. This
  changes `fit$measure` and the axis label, not any number.
* Effect modifiers coded with `contr.sum`, `contr.helmert`, `contr.poly` or a user-supplied contrast
  matrix are now covered by tests, and the coding is read from `fit$contrasts` where the fit records
  it, so a coding passed through `contrasts.arg` is recovered.
* Removed a `tests/testthat/Rplots.pdf` artefact that had been committed.

Recorded as future work rather than done: a simultaneous confidence band, a bootstrap distribution
for a selected reference location, and a continuous rather than grid-based extremum.

# svyrcs 0.5.0

## Results change when data are missing

**Knots, the reference value and the plotting range are now derived from the analytic sample** — the
rows the model is actually fitted on — rather than from the design as supplied. Previously the
spline was parameterised on one population and estimated on another whenever any model variable had
missing values.

The effect is not cosmetic. On the shipped data with missingness related to the exposure (10,617
rows to 7,748), the outer knot was 40.67 where the analytic sample gives 31.53, nine BMI units out,
and the weighted median reference was 27.35 against 25.80. The extrapolation warning added in 0.4.1
was calibrated against the full design too, so it was checking against a range the fit never saw.

Anyone comparing output across versions on a dataset with missing values should expect the numbers
to move. Point estimates and standard errors for a *given* set of knots are unchanged: fitting on
the analytic design rather than letting the model drop rows leaves them identical to 1e-15, because
survey's subsetting is proper domain estimation. What moves is the knot placement and the anchor.
Degrees of freedom fall only when a whole primary sampling unit contributes no complete case, which
is the honest answer in that situation.

* The complete-case mask covers **every** variable in the formula, not just the exposure.
* Under multiple imputation all imputations share one inclusion mask, the intersection of their
  complete cases, so that knots and coefficients mean the same thing in each. The number of rows
  this costs is reported by `print()`.
* `svyrcs::rcs(x, 4)` and `svyrcs:::rcs(x, 4)` are recognised. Previously they raised "formula must
  contain an rcs() term", which was untrue.
* A transformed exposure such as `rcs(log(bmi), 4)` is refused up front, with the workaround —
  `update(design, log_bmi = log(bmi))` — in the message.

# svyrcs 0.4.1

Findings from an adversarial verification sweep. None of them produced a wrong number, but several
let a doubtful situation pass without saying so.

* **Ordered factors work as effect modifiers.** They use polynomial contrasts, so their model matrix
  columns are `agegrp.L` and `agegrp.Q` and a level name never appears in a coefficient name;
  constructing columns from level names could never match. Group curves are now built from the
  contrast matrix, which also covers Helmert and any user-supplied contrasts. Treatment contrasts
  give bit-identical results to before.
* **A rank-deficient design now warns on every entry point instead of erroring on one.** Previously
  `svyrcs()` refused a fit that `svyrcs_curve()` would happily produce a confidence band for -- the
  same data giving opposite answers depending on how you came in. Both now warn that the design
  cannot identify the model; tests report `NA` rather than aborting the analysis.
* **Evaluating outside the observed exposure range warns.** A restricted cubic spline is linear
  beyond the outer knots, so `predict(fit, x = -50)` on a BMI model returned a number silently.
  `ref` outside the range is still an error, because anchoring there makes *every* estimate an
  extrapolation rather than just the point asked for.
* **`ref = "min"` and `ref = "max"` say when the extremum is at the range boundary.** On a monotone
  curve the "minimum-risk point" is just the end of the plotting grid, which is a property of
  `range` rather than a turning point in the data.
* A log link on a gaussian family reports `"Ratio"`, not `"RR"`: it is a ratio of means, not a risk
  ratio. Count and binary families are unaffected.
* Non-finite estimates, which perfect separation produces, are now flagged instead of landing
  silently in the curve.
* Two error messages written for specific situations were unreachable and have been fixed: colliding
  weighted quantiles now report the collision and the number of distinct values, rather than the
  number of knots that survived deduplication; and imputations with mismatched coefficients name the
  imputation at fault, rather than failing inside `vapply()`.

# svyrcs 0.4.0

* **`ggplot2` moved from `Imports` to `Suggests`.** The hard dependency chain drops from 27 packages
  to 11, six of which needed compilation. Plotting is not lost: `plot()` uses ggplot2 when it is
  installed and draws an equivalent base graphics plot when it is not.
* `plot()` gains a `backend` argument (`"auto"`, `"ggplot2"`, `"base"`) so either path can be
  requested explicitly.
* **Bug fix:** `plot()` returned the fit rather than the `ggplot`, so the documented idiom
  `plot(fit) + ggplot2::labs(...)` silently evaluated to `NULL` instead of modifying the plot.
  ggplot2's `+` accepts a non-ggplot left operand without erroring, so neither users nor `R CMD
  check` would have noticed. `plot()` now returns the `ggplot` invisibly on the ggplot2 path.
* **Breaking:** `autoplot(fit)` now needs `library(ggplot2)` first. With ggplot2 only suggested, the
  package can no longer re-export the generic; the method is registered on `ggplot2::autoplot`
  through delayed S3 registration instead. `plot(fit)` is unaffected.
* Faceted base graphics panels share one y axis. Per-panel axes would make a smaller subgroup effect
  look identical to a larger one, which is the opposite of what panelling subgroups is for.

# svyrcs 0.3.0

* `svyrcs()` accepts a multiply imputed design, built by passing a `mitools::imputationList()` to
  `survey::svydesign()`. The model is fitted in every imputation and combined by Rubin's rules.
* Knots and the reference value are computed once, as the average of the survey-weighted quantiles
  across imputations, and then held fixed. Without that the per-imputation coefficient vectors would
  not estimate the same parameters.
* **Confidence intervals use Barnard-Rubin degrees of freedom, not `mitools::MIcombine()`'s.**
  `MIcombine()` omits the small-sample correction and its degrees of freedom are unbounded — against
  a design with 31 degrees of freedom it returned values in the thousands, and `Inf`. Using those for
  a *t* quantile gives intervals far too narrow. The corrected value is bounded by the complete-data
  degrees of freedom and equals it exactly when nothing is imputed.
* The degrees of freedom are computed **per curve point**, because the between-imputation variance
  of the contrast varies along the curve.
* Joint tests get the same treatment: the Wald statistic on the pooled covariance, referred to an
  *F* distribution whose denominator applies Barnard-Rubin to the block. The Li-Raghunathan-Rubin
  denominator is not used, because it ignores the survey design.
* The curve and `predict()` gain `df` and `fmi` columns, and `print()` reports `m` and the fraction
  of missing information.
* Effect modification (`rcs(x, k) * sex`) and all four model families work with imputed designs.
* A fit with nothing to impute reproduces the complete-data fit exactly — estimates, intervals and
  tests — which is asserted in the test suite rather than merely intended.
* `mitools` and `mice` are `Suggests` only; the pooling is implemented in the package, because the
  degrees of freedom need the within- and between-imputation covariance separately and
  `MIcombine()` returns only their combination.

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
