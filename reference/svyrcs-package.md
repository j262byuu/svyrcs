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
[`rcs()`](https://j262byuu.github.io/svyrcs/reference/rcs.md) term:

    design <- svydesign(id = ~psu, strata = ~strata, weights = ~weight,
                        nest = TRUE, data = nhanes_bmi)
    fit <- svyrcs(Surv(time, event) ~ rcs(bmi, 4) + age + sex, design = design)
    summary(fit)
    plot(fit)

## Note on masking

`svyrcs` exports a function called
[`rcs()`](https://j262byuu.github.io/svyrcs/reference/rcs.md), which
masks [`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html) when
both packages are attached. The two produce numerically identical bases
(agreement to about 1e-14), so the masking does not change any result.
It does add knot storage, which is what makes
[`predict()`](https://rdrr.io/r/stats/predict.html) on new data safe.
See [`rcs()`](https://j262byuu.github.io/svyrcs/reference/rcs.md).

## See also

Useful links:

- <https://github.com/j262byuu/svyrcs>

- Report bugs at <https://github.com/j262byuu/svyrcs/issues>

## Author

**Maintainer**: Xiaoyu Zong <j262byuu@gmail.com>

Authors:

- Xiaoyu Zong <j262byuu@gmail.com>
