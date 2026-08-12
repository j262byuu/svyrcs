# svyrcs: Restricted Cubic Splines for Complex Survey Data

Fits and summarises restricted cubic spline (RCS) models for complex
survey data such as the National Health and Nutrition Examination Survey
(NHANES). Splines are combined with survey-weighted Cox and generalised
linear models from the 'survey' package, and exposure-response curves
are estimated as contrasts against a reference value with design-based
variance – Taylor linearisation for ordinary designs, replicate-weight
variance for replicate designs – and confidence intervals on the survey
degrees of freedom by default. Provides tests of overall association and
of non-linearity, several ways of choosing the reference value, and
publication-ready plots.

## Getting started

Build a survey design with
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html),
then call
[`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md) with
a model formula containing one
[`rcspline()`](https://j262byuu.github.io/svyrcs/reference/rcspline.md)
term:

    design <- svydesign(id = ~psu, strata = ~strata, weights = ~weight,
                        nest = TRUE, data = nhanes_bmi)
    fit <- svyrcs(Surv(time, event) ~ rcspline(bmi, 4) + age + sex, design = design)
    summary(fit)
    plot(fit)

## Using this package alongside `rms`

The spline function here is
[`rcspline()`](https://j262byuu.github.io/svyrcs/reference/rcspline.md).
It was called `rcs()` before 0.7.0, which masked
[`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html) whenever both
packages were attached; the name was changed because that masking could
not be made safe. A `cph()`, `lrm()` or `ols()` model written with a
bare `rcs()` would silently pick up this basis, fit correctly, and then
lose the `Nonlinear` row from
[`anova()`](https://rdrr.io/r/stats/anova.html) – with no error and no
warning, and that row is usually the reason for fitting a spline at all.

Since the rename there is no collision, and no attach order to remember.
Use [`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html) in `rms`
models and
[`rcspline()`](https://j262byuu.github.io/svyrcs/reference/rcspline.md)
here. Note that the converse is not symmetric: passing
[`rcspline()`](https://j262byuu.github.io/svyrcs/reference/rcspline.md)
into an `rms` model still fits and still loses the `Nonlinear` test,
because the basis carries this package's knot attributes rather than
`rms`'s design attributes.

## See also

Useful links:

- <https://github.com/j262byuu/svyrcs>

- Report bugs at <https://github.com/j262byuu/svyrcs/issues>

## Author

**Maintainer**: Xiaoyu Zong <j262byuu@gmail.com>

Authors:

- Xiaoyu Zong <j262byuu@gmail.com>
