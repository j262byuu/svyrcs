# svyrcs

Restricted cubic splines for complex survey data such as NHANES, with
design-based variance estimation.

Splines and survey designs do not normally meet. `rms` gives you
restricted cubic splines but knows nothing about sampling weights,
clustering or stratification; `survey` handles all three but has no
spline helper. `svyrcs` fits the spline *inside* a `survey` model, so
point estimates use the weights and the confidence band uses Taylor
linearisation on the design degrees of freedom.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("j262byuu/svyrcs")
```

## Quick start

``` r

library(svyrcs)
library(survey)
library(survival)

design <- svydesign(
  id = ~psu, strata = ~strata, weights = ~weight,
  nest = TRUE, data = nhanes_bmi
)

fit <- svyrcs(Surv(time, event) ~ rcs(bmi, 4) + age + sex + race, design = design)
fit
#> Restricted cubic spline on a complex survey design
#>
#>   Model      survey-weighted Cox proportional hazards
#>   Outcome    Surv(time, event)
#>   Exposure   bmi, 4 knots at 19.78, 25.22, 29.71, 40.67 (weighted quantiles)
#>   Sample     10,617 observations, 1,893 events, 31 design df
#>   Reference  bmi = 27.35 (weighted median), HR = 1
#>
#>   Overall association  F = 43.58 on 3 and 31 df,  p = 3.12e-11
#>   Non-linearity        F = 51.70 on 2 and 31 df,  p = 1.34e-10

plot(fit)                          # ggplot curve with confidence band
summary(fit)                       # turning points, where the band excludes the null
predict(fit, x = c(25, 30, 35))    # estimates at values you choose
```

Anchor the curve wherever the analysis calls for it:

``` r

svyrcs(Surv(time, event) ~ rcs(bmi, 4) + age, design = design, ref = "min")
```

Estimate one curve per subgroup, from a single model, with a test of
whether they differ:

``` r

by_sex <- svyrcs(Surv(time, event) ~ rcs(bmi, 4) * sex + age, design = design)
by_sex$tests$interaction   # p for interaction
by_sex$tests$shape         # does the *curvature* differ, or is the curve just shifted?
plot(by_sex)               # one coloured curve per group; facet = TRUE to panel them
```

Any outcome `survey` can fit:

| outcome      | call                                        | reports         |
|--------------|---------------------------------------------|-----------------|
| survival     | `Surv(time, event) ~ rcs(x, 4)`             | hazard ratio    |
| binary       | `y ~ rcs(x, 4)`, `family = quasibinomial()` | odds ratio      |
| rate / count | `y ~ rcs(x, 4)`, `family = quasipoisson()`  | rate ratio      |
| continuous   | `y ~ rcs(x, 4)`                             | mean difference |

## What this does that a hand-rolled script does not

This package grew out of a tutorial script. Turning it into a package
fixed several things that are easy to get wrong when doing this by hand:

| the usual approach | what `svyrcs` does | why it matters |
|----|----|----|
| build a prediction grid with every covariate fixed at a reference value, then [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) | contrast the spline columns only | covariates cancel algebraically, so the grid is unnecessary; it also breaks when a covariate is a factor with a level missing from the grid |
| find spline coefficients by `grep`-ing the variable name | exact match on the term label | a covariate called `bmi_group` would otherwise be swept into the estimate for `bmi` |
| chi-square *p*-values | design-based *F* on `(df1, degf)` | with 31 design degrees of freedom the chi-square test is anti-conservative; the *F* test matches [`survey::regTermTest()`](https://rdrr.io/pkg/survey/man/regTermTest.html) |
| normal quantile for the confidence band | *t* quantile on the design df | a normal quantile gives intervals that are too narrow |
| knots and reference at unweighted quantiles | survey-weighted quantiles by default | in a design-based analysis the population median is the weighted one; `weighted_knots = FALSE` restores the old behaviour |
| [`rms::rcs()`](https://rdrr.io/pkg/rms/man/rms.trans.html) plus the `datadist` global option | knots stored on the fitted model | [`predict()`](https://rdrr.io/r/stats/predict.html) on new data can never silently re-derive different knots, and there is no global state to keep in sync |
| subgroup curves from separate stratified fits | one interaction model, curves via a selection matrix | separate fits cannot produce a *p* for interaction, and estimate the covariate effects independently; the matrix form also keeps `Cov(main, interaction)` in the standard error, which adding the two variances would drop |

## Documentation

[`vignette("svyrcs")`](https://j262byuu.github.io/svyrcs/articles/svyrcs.md)
walks through the whole workflow on the shipped NHANES data.

`nhanes_bmi` is adults from NHANES 2005–2008 linked to the NCHS 2019
Linked Mortality File, both cycles complete so the design is valid. It
is public domain and ships with the package, so every example runs on
real survey data.

## Citation

If this package was useful, please cite it:

``` r

citation("svyrcs")
```

## Contributing

Issues and pull requests are welcome at
<https://github.com/j262byuu/svyrcs>. If you think a number is wrong, an
issue with a reproducible example is the most useful thing you can send.

## References

1.  Harrell FE Jr. *Regression Modeling Strategies*. 2nd
    ed. Springer; 2015.
    [doi:10.1007/978-3-319-19425-7](https://doi.org/10.1007/978-3-319-19425-7)
2.  Lumley T. *Complex Surveys: A Guide to Analysis Using R*.
    Wiley; 2010.
    [doi:10.1002/9780470580066](https://doi.org/10.1002/9780470580066)
3.  Binder DA. On the variances of asymptotically normal estimators from
    complex surveys. *International Statistical Review*.
    1983;51(3):279–292.
    [doi:10.2307/1402588](https://doi.org/10.2307/1402588)
4.  CDC. [NHANES
    tutorials](https://wwwn.cdc.gov/nchs/nhanes/tutorials/default.aspx)

## License

MIT © svyrcs authors
