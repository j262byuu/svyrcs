# Package index

## Fitting

Fit a spline model to a survey design and get the curve, the tests and
the reference.

- [`svyrcs()`](https://j262byuu.github.io/svyrcs/reference/svyrcs.md) :
  Restricted cubic spline analysis of complex survey data
- [`svyrcs_curve()`](https://j262byuu.github.io/svyrcs/reference/svyrcs_curve.md)
  : Exposure-response curve from a fitted survey model

## The spline basis

Building and inspecting the restricted cubic spline basis.

- [`rcs()`](https://j262byuu.github.io/svyrcs/reference/rcs.md) :
  Restricted cubic spline basis
- [`rcs_knots()`](https://j262byuu.github.io/svyrcs/reference/rcs_knots.md)
  : Knot locations for a restricted cubic spline
- [`makepredictcall(`*`<svyrcs_basis>`*`)`](https://j262byuu.github.io/svyrcs/reference/makepredictcall.svyrcs_basis.md)
  : Reuse fitted knots when predicting on new data

## Results

Looking at, summarising and plotting a fit.

- [`summary(`*`<svyrcs>`*`)`](https://j262byuu.github.io/svyrcs/reference/summary.svyrcs.md)
  : Summarise a survey restricted cubic spline fit
- [`autoplot(`*`<svyrcs>`*`)`](https://j262byuu.github.io/svyrcs/reference/plot.svyrcs.md)
  [`plot(`*`<svyrcs>`*`)`](https://j262byuu.github.io/svyrcs/reference/plot.svyrcs.md)
  : Plot a survey restricted cubic spline curve
- [`predict(`*`<svyrcs>`*`)`](https://j262byuu.github.io/svyrcs/reference/predict.svyrcs.md)
  : Effect estimates at specific exposure values
- [`references`](https://j262byuu.github.io/svyrcs/reference/references.md)
  : Reference values for an exposure-response curve

## Data

- [`nhanes_bmi`](https://j262byuu.github.io/svyrcs/reference/nhanes_bmi.md)
  : NHANES 2005-2008 with linked mortality follow-up

## Package

- [`svyrcs-package`](https://j262byuu.github.io/svyrcs/reference/svyrcs-package.md)
  : svyrcs: Restricted Cubic Splines for Complex Survey Data
