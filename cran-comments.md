# CRAN submission comments

## Submission

svyrcs 0.4.0 — new submission.

The package fits restricted cubic splines inside `survey` models, so that exposure-response curves
for complex survey data (NHANES and similar) are estimated with sampling weights and design-based
variance. It fills a gap between `rms`, which provides the spline basis but does not handle survey
designs, and `survey`, which handles the design but has no spline helper.

## Test environments

Checked with `--as-cran` on:

| Platform | R |
|---|---|
| macOS 15, aarch64-apple-darwin23 (GitHub Actions) | 4.6.1 |
| Windows Server 2022, x86_64-w64-mingw32 (GitHub Actions) | 4.6.1 |
| Ubuntu 24.04, x86_64-pc-linux-gnu (GitHub Actions) | 4.6.1 |
| Ubuntu 24.04, x86_64-pc-linux-gnu (GitHub Actions) | R-devel (2026-06-21 r90185) |
| Ubuntu 24.04, x86_64-pc-linux-gnu (GitHub Actions) | 4.5.3 |
| Windows 11, x86_64-w64-mingw32 (local) | 4.5.2 |

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Xiaoyu Zong <j262byuu@gmail.com>'
  New submission

All five GitHub Actions platforms reported `Status: OK` with no notes at all. The local run reports
one additional note, `Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being
installed`, which is an artefact of that machine: both files were verified separately with pandoc
and convert cleanly.

The local run also uses `--no-manual`, since no LaTeX is installed there. The manual builds on the
GitHub Actions platforms.

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

Examples run in well under 5 seconds each. Three survival examples, which fit survey-weighted Cox
models, are wrapped in `\donttest{}` to keep the total example time low; they run cleanly under
`--run-donttest`, which was included in every check above.

The test suite is 634 tests and takes about 60 seconds. Tests that need a `Suggests` package are
guarded with `skip_if_not_installed()`.

## Dependencies

`Imports` is limited to `rlang`, `stats` and `survey` — 11 packages including their recursive
dependencies. `ggplot2`, `survival`, `rms`, `mitools`, `mice`, `haven`, `knitr`, `rmarkdown` and
`testthat` are in `Suggests` and are used only where guarded.

`ggplot2` is suggested rather than imported even though plotting is a headline feature: `plot()`
falls back to an equivalent base graphics implementation when ggplot2 is absent, so nothing is lost.
The `autoplot()` method is attached to `ggplot2::autoplot` by delayed S3 registration
(`S3method(ggplot2::autoplot, svyrcs)`), so the generic is never required at load time.

The package deliberately does not depend on `rms`, even though it implements the same spline basis:
the two agree to about 1e-14, and avoiding the dependency also avoids `rms`'s `datadist` global
option.

## Note on masking

The package exports a function named `rcs()`, which masks `rms::rcs()` when both packages are
attached. This is intentional and documented in `?rcs` and in the package-level help. The two produce
numerically identical bases, so the masking cannot change a result; `svyrcs::rcs()` additionally
stores its knots on the returned object and registers a `makepredictcall()` method, which is what
makes prediction on new data safe.

## Downstream dependencies

None — this is a new package.
