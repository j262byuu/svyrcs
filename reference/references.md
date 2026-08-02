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

  The exposure value at which the fitted curve is lowest: the
  minimum-risk point.

- `"max"`:

  The exposure value at which the fitted curve is highest.

Weighted quantiles and means are computed from the survey design with
[`survey::svyquantile()`](https://rdrr.io/pkg/survey/man/svyquantile.html)
and
[`survey::svymean()`](https://rdrr.io/pkg/survey/man/surveysummary.html),
so the reference is a population quantity rather than a property of the
unweighted sample. `"min"` and `"max"` are found by evaluating the curve
on a dense grid against a provisional (weighted median) reference;
because the reference only shifts the curve, the located point does not
depend on that provisional choice.
