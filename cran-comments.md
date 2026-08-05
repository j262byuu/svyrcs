# CRAN submission comments

## Submission

svyrcs 0.6.4 — new submission.

The package fits restricted cubic splines inside `survey` models, so that exposure-response curves
for complex survey data (NHANES and similar) are estimated with sampling weights and design-based
variance. It fills a gap between `rms`, which provides the spline basis but does not handle survey
designs, and `survey`, which handles the design but has no spline helper.

## Test environments

Checked with `--as-cran` on:

| Platform | R |
|---|---|
| macOS Tahoe 26.5.2, aarch64-apple-darwin23 (GitHub Actions) | 4.6.1 (2026-06-24) |
| Windows Server 2025, x86_64-w64-mingw32 (GitHub Actions) | 4.6.1 (2026-06-24 ucrt) |
| Ubuntu 24.04.4 LTS, x86_64-pc-linux-gnu (GitHub Actions) | 4.6.1 (2026-06-24) |
| Ubuntu 24.04.4 LTS, x86_64-pc-linux-gnu (GitHub Actions) | R-devel (2026-06-21 r90185) |
| Ubuntu 24.04.4 LTS, x86_64-pc-linux-gnu (GitHub Actions) | 4.5.3 (2026-03-11) |
| Ubuntu 24.04.4 LTS, x86_64-pc-linux-gnu, hard dependencies only (GitHub Actions) | 4.6.1 (2026-06-24) |
| Windows 11, x86_64-w64-mingw32 (local) | 4.5.2 |

All six GitHub Actions jobs reported `Status: OK`.

The hard-dependencies-only run installs nothing from `Suggests` beyond what the test suite itself
needs, and asserts that ggplot2 is absent before checking, so the base graphics fallback is
exercised rather than assumed. It reported `Status: OK` with 849 tests passing and 19 skipped; the
other five jobs ran the full suite of 1037 tests with nothing skipped.

One job additionally installs TinyTeX and HTML Tidy and drops `--no-manual`, so both the PDF and the
HTML manual are built and checked. A post-check step asserts that those two lines are present and
passing, since a skipped check and a passing check are otherwise indistinguishable from the job's
exit status.

## R CMD check results

0 errors | 0 warnings | 1 note on the CI platforms:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Xiaoyu Zong <j262byuu@gmail.com>'
  New submission

The local run reports one further note, `checking for future file timestamps ... unable to verify
current time`, which appears when the machine cannot reach a time server. It is transient and
environmental.

## Shipped data

`nhanes_bmi` (134 KB, xz-compressed) is an extract of the National Health and Nutrition Examination
Survey 2005-2006 and 2007-2008 public use files, joined to the National Center for Health Statistics
2019 Linked Mortality File. Both are works of the United States Government and are in the public
domain, so redistribution is permitted. The script that builds the extract is included in
`data-raw/`, and the sources are documented in `?nhanes_bmi`.

Two complete survey cycles are shipped rather than a sample of rows: sampling rows or primary
sampling units would invalidate the design, and every example and test in the package relies on the
design being valid.

## Examples and tests

Each example runs in under 5 seconds. Three survival examples, which fit survey-weighted Cox models,
are wrapped in `\donttest{}`, as are the secondary plotting calls in `?plot.svyrcs` — the base
graphics fallback and the ggplot2 modification. All of them run cleanly under `--run-donttest`, which
was included in every check above.

The test suite is 1037 tests and takes about three minutes. Tests that need a `Suggests` package are
guarded with `skip_if_not_installed()`.

## Dependencies

`Imports` is limited to `rlang`, `stats` and `survey (>= 4.0)` — 11 packages to install, counting
their recursive dependencies and excluding the base packages that ship with R (18 including them).

`ggplot2`, `haven`, `Hmisc`, `knitr`, `mice`, `mitml`, `mitools`,
`rmarkdown`, `rms`, `survival` and `testthat` are in `Suggests` and are used only where guarded.

`ggplot2` is suggested rather than imported even though plotting is a headline feature: `plot()`
falls back to an equivalent base graphics implementation when ggplot2 is absent, so nothing is lost.
The `autoplot()` method is attached to `ggplot2::autoplot` by delayed S3 registration
(`S3method(ggplot2::autoplot, svyrcs)`), so the generic is never required at load time.

A dedicated continuous integration job installs only the hard dependencies plus what the test suite
itself needs, and asserts that ggplot2 really is absent before running, so the base graphics fallback
is exercised in a real environment rather than only under mocking.

`Hmisc` and `rms` are suggested purely so that the test suite can check this package's spline basis
and knot placement against Harrell's reference implementation; neither is used at run time. The
package deliberately does not depend on `rms`, even though it implements the same spline basis: the
two agree to about 1e-14, and avoiding the dependency also avoids `rms`'s `datadist` global option.

## Note on masking

The package exports a function named `rcs()`, which masks `rms::rcs()` when both packages are
attached. This is intentional and documented in `?rcs` and in the package-level help. The two produce
numerically identical bases given the same knots, and `rcs_knots()` reproduces
`Hmisc::rcspline.eval()`'s default placement exactly, so the masking cannot change a result.
`svyrcs::rcs()` additionally stores its knots on the returned object and registers a
`makepredictcall()` method, which is what makes prediction on new data safe.

## Downstream dependencies

None — this is a new package.
