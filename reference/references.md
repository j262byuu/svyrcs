# Reference values for an exposure-response curve

A restricted cubic spline curve is only identified up to a reference
value: every estimate is a contrast against some exposure level, which
is where the curve passes through the null. `svyrcs` supports the
reference choices commonly seen in the survey literature.

## Available methods

- a number:

  Used as given. Must lie inside the observed exposure range.

- `"median"`:

  The survey-weighted median of the exposure. The default.

- `"mean"`:

  The survey-weighted mean.

- `"quantile"`:

  The survey-weighted quantile at `ref_prob`.

- `"min"`:

  The grid point at which the fitted curve is lowest.

- `"max"`:

  The grid point at which the fitted curve is highest.

Weighted quantiles and means are computed from the survey design with
[`survey::svyquantile()`](https://rdrr.io/pkg/survey/man/svyquantile.html)
and
[`survey::svymean()`](https://rdrr.io/pkg/survey/man/surveysummary.html),
so the reference is a population quantity rather than a property of the
unweighted sample. `"min"` and `"max"` are found by evaluating the curve
on a dense grid against a provisional (weighted median) reference;
because the reference only shifts the curve, the located point does not
depend on that provisional choice.

## Selected references are descriptive

`"min"` and `"max"` choose an anchor **from the same data that produced
the curve**, and that choice carries sampling uncertainty which the
confidence band does not propagate. Two consequences worth stating
plainly:

- At the selected point the estimate equals the null with a standard
  error of zero. That is an artefact of anchoring the contrast there,
  not evidence that the null holds at that exposure.

- The location itself is not estimated with the precision the plot
  suggests. In another sample the minimum could sit some distance away,
  particularly where the curve is flat near its bottom.

Treat them as a presentational choice — "excess risk relative to the
healthiest observed value" — rather than as an estimate of an optimum.
Reporting a confidence interval for the location would need a bootstrap
or replicate-weight distribution over refits, which this package does
not currently provide.

Both are also grid searches, so the answer can shift slightly with `n`.
