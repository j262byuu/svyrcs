# NHANES 2005-2008 with linked mortality follow-up

Adults examined in the 2005-2006 and 2007-2008 cycles of the National
Health and Nutrition Examination Survey, linked to the National Center
for Health Statistics 2019 Linked Mortality File. Supplied so that every
example and test in this package runs on a real complex survey rather
than on simulated data.

## Usage

``` r
nhanes_bmi
```

## Format

A data frame with 10,617 rows and 16 columns:

- seqn:

  NHANES respondent sequence number.

- cycle:

  Survey cycle, `"2005-2006"` or `"2007-2008"`.

- psu:

  Masked variance unit (`SDMVPSU`), the primary sampling unit.

- strata:

  Masked variance stratum (`SDMVSTRA`).

- weight:

  MEC examination weight for the pooled two cycles, i.e. `WTMEC2YR / 2`.

- age:

  Age in years at screening (top-coded at 80).

- sex:

  `"Male"` or `"Female"`.

- race:

  Race and Hispanic origin, four levels.

- bmi:

  Body mass index, kg/m2.

- tchol:

  Total cholesterol, mg/dL.

- hdl:

  HDL cholesterol, mg/dL.

- glucose:

  Fasting plasma glucose, mg/dL. Only measured in the fasting subsample,
  so mostly missing; a fasting-subsample analysis needs its own weight,
  which is not included here.

- statin_use:

  1 if a statin was reported in the prescription medication interview.

- high_chol:

  1 if `tchol` is 240 mg/dL or more.

- time:

  Follow-up from the examination to death or 31 December 2019, in years.

- event:

  1 if the participant died during follow-up, 0 if censored.

## Source

NHANES 2005-2006 and 2007-2008 public use files, and the NCHS 2019
Linked Mortality File, <https://www.cdc.gov/nchs/nhanes/>. Both are in
the public domain as works of the US Government. The script that builds
this extract is `data-raw/make_nhanes_bmi.R` in the package sources.

## Details

Both survey cycles are included in full, so all 31 sampling strata and
all 62 primary sampling units are intact and the design has 31 degrees
of freedom. Sampling rows or PSUs would have made the design invalid, so
nothing was sampled: the reduction in size comes entirely from keeping
two cycles instead of nine.

Adults who were examined and not pregnant, eligible for mortality
linkage, and with non-missing body mass index, age, sex, race and
follow-up time. Analyses of the fasting subsample (`glucose`) would need
the fasting weight and are not supported by the weight shipped here.

## Examples

``` r
design <- survey::svydesign(
  id = ~psu, strata = ~strata, weights = ~weight,
  nest = TRUE, data = nhanes_bmi
)
survey::degf(design)
#> [1] 31
summary(nhanes_bmi$bmi)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   13.36   24.00   27.62   28.66   31.93  130.21 
```
