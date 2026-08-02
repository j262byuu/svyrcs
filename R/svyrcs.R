## Walk a formula's language tree and rewrite the rcs() call so it carries `knots` literally.
##
## The knots have to be inlined rather than referred to through a symbol: survey::svycoxph() and
## survey::svyglm() both reject a formula whose all.vars() are not all columns of the design. Inlining
## also makes the fitted model self-contained -- predict(), regTermTest() and saveRDS() all work
## without any environment tagging along.
rewrite_rcs_call <- function(e, knots) {
  if (is.call(e)) {
    if (identical(as.character(e[[1L]])[1L], "rcs")) {
      out <- e[1L:2L]
      out$knots <- knots
      return(out)
    }
    for (i in seq_along(e)[-1L]) e[[i]] <- rewrite_rcs_call(e[[i]], knots)
  }
  e
}

## Shorten knots for the sake of readable coefficient names. Knot placement is a modelling choice,
## not an estimate, so rounding it is harmless -- as long as the rounded values are then used
## everywhere, which they are: these are the knots stored on the returned object.
tidy_knots <- function(knots, digits = 7L) {
  out <- signif(knots, digits)
  if (anyDuplicated(out) || length(unique(out)) < length(knots)) knots else out
}

## Every rcs() call appearing anywhere in a formula.
collect_rcs_calls <- function(e, acc = list()) {
  if (is.call(e)) {
    if (identical(as.character(e[[1L]])[1L], "rcs")) return(c(acc, list(e)))
    for (i in seq_along(e)[-1L]) acc <- collect_rcs_calls(e[[i]], acc)
  }
  acc
}

## Build the `tests` component.
##
## Without a modifier this is exactly the 0.1.0 pair. With one, `overall` and `nonlinear` widen to
## the joint test over the spline block *and* its interactions -- "is the exposure associated with
## the outcome in any group?" -- and the per-group and interaction tests are added alongside. Keeping
## the same two names in both shapes means code that reads fit$tests$overall keeps working.
assemble_tests <- function(beta, V, term, modifier, degf) {
  spline_cols <- spline_colnames(term$label, term$nk)

  if (is.null(modifier)) {
    L <- selection_matrix(names(beta), spline_cols, NULL)
    return(group_tests(beta, V, L, degf))
  }

  int_cols <- unlist(modifier$columns, use.names = FALSE)
  joint <- c(spline_cols, int_cols)
  joint_nl <- c(spline_cols[-1L],
                unlist(lapply(modifier$columns, function(cols) cols[-1L]), use.names = FALSE))

  by_group <- lapply(modifier$levels, function(g) {
    group_tests(beta, V, selection_matrix(names(beta), spline_cols, modifier$columns[[g]]), degf)
  })
  names(by_group) <- modifier$levels

  itests <- interaction_tests(beta, V, modifier, degf)

  list(
    overall = wald_block(beta[joint], V[joint, joint, drop = FALSE], degf, "overall association"),
    nonlinear = wald_block(beta[joint_nl], V[joint_nl, joint_nl, drop = FALSE], degf,
                           "non-linearity"),
    interaction = itests$interaction,
    shape = itests$shape,
    by_group = by_group
  )
}

## Detect a Surv() response, whether written bare or namespace-qualified as survival::Surv().
is_surv_call <- function(e) {
  if (!is.call(e)) return(FALSE)
  fn <- e[[1L]]
  nm <- if (is.call(fn) && identical(as.character(fn[[1L]])[1L], "::")) {
    as.character(fn[[3L]])
  } else {
    as.character(fn)[1L]
  }
  identical(nm, "Surv")
}

#' Restricted cubic spline analysis of complex survey data
#'
#' Fits a restricted cubic spline model to a complex survey design and returns the
#' exposure-response curve, its confidence band, tests of overall association and of non-linearity,
#' and the reference value the curve is anchored to.
#'
#' @param formula A model formula containing exactly one [rcs()] term for the exposure, plus any
#'   covariates. A `Surv()` response selects a survey-weighted Cox model; anything else is fitted
#'   with [survey::svyglm()]. Crossing the spline with a grouping variable, `rcs(x, 4) * sex`,
#'   estimates one curve per group and tests whether they differ; see the section below.
#' @param design A survey design object, as built by [survey::svydesign()]. Subset it with
#'   [subset()] before passing it in; the degrees of freedom are taken from the design you pass.
#' @param family Family for the [survey::svyglm()] path. Use `quasibinomial()` for a binary outcome
#'   (odds ratios), `quasipoisson()` for counts or rates (rate ratios), or the default
#'   `gaussian()` for a continuous outcome (mean differences). Ignored for Cox models.
#' @param ref Reference value the curve is anchored to: a number, or one of `"median"` (default),
#'   `"mean"`, `"quantile"`, `"min"` or `"max"`. See [references].
#' @param ref_prob Probability used when `ref = "quantile"`.
#' @param weighted_knots If `TRUE` (default), a knot *count* is turned into knot locations using
#'   survey-weighted quantiles of the exposure, so knot placement reflects the population rather than
#'   the unweighted sample. Set to `FALSE` for unweighted quantiles. Explicit knot locations are
#'   always used as given.
#' @param n Number of points on the estimated curve.
#' @param range Length-2 numeric giving the exposure range to estimate over. Defaults to the 1st and
#'   99th survey-weighted percentiles.
#' @param level Confidence level for the band.
#' @param ... Passed on to [survey::svycoxph()] or [survey::svyglm()].
#'
#' @return An object of class `svyrcs`, a list with components:
#'   \describe{
#'     \item{curve}{Data frame of `x`, `estimate`, `conf.low`, `conf.high`, `se`.}
#'     \item{ref}{The reference value and how it was chosen.}
#'     \item{knots}{Knot locations actually used.}
#'     \item{tests}{Overall and non-linearity Wald tests, as *F* and chi-square.}
#'     \item{measure}{`"HR"`, `"OR"`, `"RR"` or `"Difference"`.}
#'     \item{model}{The fitted `svycoxph`/`svyglm` object, for your own diagnostics.}
#'     \item{degf, n, nevents, level, var, call}{Analysis metadata.}
#'   }
#'
#' @section Subgroups and effect modification:
#' Writing the exposure crossed with a grouping variable, `rcs(bmi, 4) * sex`, fits **one** model and
#' returns one curve per level of the modifier. Because it is a single model, the covariate effects
#' are shared across groups and a genuine interaction test is available; fitting each subgroup
#' separately gives neither. `sex * rcs(bmi, 4)` is equivalent.
#'
#' The modifier must be a factor, character or logical variable. Two extra tests are then reported:
#' \describe{
#'   \item{Interaction}{Are the spline-by-group terms jointly zero, i.e. does the association differ
#'     between groups at all?}
#'   \item{Shape interaction}{The same test without the linear term: does the *curvature* differ, as
#'     opposed to the whole curve being shifted?}
#' }
#' plus the overall and non-linearity tests within each group.
#'
#' All groups share one reference value, so that the curves are comparable; `ref = "min"` and
#' `ref = "max"` are located on the reference level's curve, which the printed output states.
#'
#' @section Why not just use `rms`:
#' `rms::rcs()` gives the same basis, but `rms` fitting functions do not know about sampling weights,
#' clustering or stratification, and `survey` has no spline helper. `svyrcs` bridges the two: the
#' spline is fitted inside a `survey` model, so point estimates use the weights and the confidence
#' band uses Taylor linearisation on the design degrees of freedom.
#'
#' @seealso [svyrcs_curve()] for curves from a model you fitted yourself, [rcs()] for the basis,
#'   [references] for reference values.
#'
#' @examples
#' design <- survey::svydesign(
#'   id = ~psu, strata = ~strata, weights = ~weight,
#'   nest = TRUE, data = nhanes_bmi
#' )
#'
#' # continuous outcome: mean difference in total cholesterol across the BMI range
#' fit <- svyrcs(tchol ~ rcs(bmi, 4) + age + sex, design = design)
#' fit
#'
#' # survival outcome: all-cause mortality, anchored at the minimum-risk BMI
#' \donttest{
#' library(survival)
#' fit_hr <- svyrcs(Surv(time, event) ~ rcs(bmi, 4) + age + sex,
#'                  design = design, ref = "min")
#' summary(fit_hr)
#' plot(fit_hr)
#'
#' # binary outcome: odds of high total cholesterol
#' fit_or <- svyrcs(high_chol ~ rcs(bmi, 4) + age + sex, design = design,
#'                  family = quasibinomial())
#'
#' # one curve per subgroup, with a test of whether they differ
#' fit_by_sex <- svyrcs(Surv(time, event) ~ rcs(bmi, 4) * sex + age, design = design)
#' fit_by_sex$tests$interaction$p_F
#' plot(fit_by_sex)
#' }
#'
#' @export
svyrcs <- function(formula, design, family = NULL, ref = "median", ref_prob = 0.5,
                   weighted_knots = TRUE, n = 200, range = NULL, level = 0.95, ...) {
  cl <- match.call()

  if (!inherits(design, c("survey.design", "svyrep.design"))) {
    stop_svyrcs("`design` must be a survey design object from survey::svydesign() or ",
                "survey::svrepdesign(); got ", class(design)[1L])
  }
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop_svyrcs("`formula` must be a two-sided model formula, such as ",
                "outcome ~ rcs(bmi, 4) + age")
  }

  rcs_calls <- collect_rcs_calls(formula[[3L]])
  if (length(rcs_calls) == 0L) {
    stop_svyrcs("`formula` must contain an rcs() term for the exposure, for example ",
                "outcome ~ rcs(bmi, 4) + age")
  }
  if (length(rcs_calls) > 1L) {
    stop_svyrcs("`formula` contains ", length(rcs_calls), " rcs() terms; svyrcs models one ",
                "exposure at a time. Fit separate models, or use svyrcs_curve() on a model you ",
                "fit yourself.")
  }

  rcs_call <- rcs_calls[[1L]]
  var <- deparse1(rcs_call[[2L]])
  if (is.null(design$variables[[var]])) {
    stop_svyrcs("the exposure '", var, "' is not a variable in `design`")
  }
  knot_arg <- if (!is.null(rcs_call$knots)) rcs_call$knots else
    if (length(rcs_call) >= 3L) rcs_call[[3L]] else 4
  knot_spec <- eval(knot_arg, envir = design$variables, enclos = environment(formula))

  ## Resolve the knots up front so they are fixed for fitting, prediction and plotting alike.
  xvals <- as.numeric(design$variables[[var]])
  knots <- if (length(knot_spec) == 1L && weighted_knots) {
    if (!is_count(knot_spec)) {
      stop_svyrcs("`knots` of length 1 is read as the number of knots and must be a whole number ",
                  "between 3 and 7; got ", format(knot_spec))
    }
    nk <- as.integer(round(knot_spec))
    probs <- harrell_knot_probs(nk)
    kn <- unname(weighted_quantile(var, design, probs))
    rcs_knots(xvals, kn, var = var)
  } else {
    rcs_knots(xvals, knot_spec, var = var)
  }
  knots <- tidy_knots(knots)

  ## Rewrite the formula so the fit uses these exact knots.
  fit_formula <- formula
  fit_formula[[3L]] <- rewrite_rcs_call(formula[[3L]], knots)

  is_cox <- is_surv_call(formula[[2L]])
  if (is_cox) {
    if (!requireNamespace("survival", quietly = TRUE)) {
      stop_svyrcs("a Surv() response needs the 'survival' package; install it with ",
                  "install.packages(\"survival\")")
    }
    if (!is.null(family)) {
      warning("`family` is ignored for a Surv() response", call. = FALSE)
    }
    model <- survey::svycoxph(fit_formula, design = design, ...)
  } else {
    model <- survey::svyglm(fit_formula, design = design,
                            family = family %||% stats::gaussian(), ...)
  }

  degf <- survey::degf(design)
  term <- find_rcs_term(model, var)
  spline_coef_index(model, term)  # validates the main effects are present
  modifier <- find_modifier(model, term)
  beta <- stats::coef(model)
  V <- stats::vcov(model)

  curve <- svyrcs_curve(model, var = var, ref = ref, ref_prob = ref_prob, n = n, range = range,
                        level = level, design = design, degf = degf)
  tests <- assemble_tests(beta, V, term, modifier, degf)
  meas <- effect_measure(model)

  ## nobs() on a coxph object returns the number of *events* (its effective sample size for AIC),
  ## not the number of observations, so ask the model directly on the Cox path.
  n_used <- if (is_cox) as.integer(model$n) else as.integer(stats::nobs(model))
  n_design <- nrow(design$variables)
  nevents <- if (is_cox) as.integer(model$nevent) else NA_integer_

  structure(
    list(
      curve = curve,
      ref = list(value = attr(curve, "ref"), method = attr(curve, "ref_method")),
      knots = knots,
      nk = length(knots),
      groups = if (is.null(modifier)) NULL else {
        list(var = modifier$var, levels = modifier$levels, ref_level = modifier$ref_level)
      },
      var = var,
      label = paste0("rcs(", var, ", ", length(knots), ")"),
      ## The term exactly as it appears in the fitted model, so that
      ## survey::regTermTest(fit$model, reformulate(fit$term_label)) just works.
      term_label = term$label,
      measure = meas$measure,
      null = meas$null,
      exponentiate = meas$exponentiate,
      tests = tests,
      degf = degf,
      n = n_used,
      n_dropped = as.integer(n_design - n_used),
      nevents = nevents,
      level = level,
      weighted_knots = isTRUE(weighted_knots) && length(knot_spec) == 1L,
      model = model,
      formula = formula,
      call = cl
    ),
    class = "svyrcs"
  )
}
