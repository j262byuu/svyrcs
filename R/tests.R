## Joint Wald test on a coefficient block.
##
## Reported both as a chi-square and as a design-based F. With a survey design the chi-square
## reference distribution ignores the fact that the variance is estimated from a modest number of
## PSUs, so it is anti-conservative; F on (df1, degf) is the survey convention and reproduces
## survey::regTermTest() exactly.
## `mi`, when supplied, must already hold the within (`Ubar`) and between (`B`) covariance blocks in
## the same basis as `V`, so that every caller -- whole spline block, per-group selection, the
## interaction columns -- transforms once and this function stays basis-agnostic.
wald_block <- function(beta, V, degf, what, mi = NULL) {
  df1 <- length(beta)
  if (df1 < 1L) {
    return(list(chisq = NA_real_, F = NA_real_, df1 = 0L, df2 = degf,
                p_chisq = NA_real_, p_F = NA_real_, fmi = if (is.null(mi)) NULL else 0))
  }
  ## Under imputation the joint test is D1, not a Wald statistic on the pooled covariance. The
  ## statistic, its denominator degrees of freedom and the reported chi-square all come from there.
  if (!is.null(mi)) {
    d1 <- d1_test(beta, mi$Ubar, mi$B, mi$m, mi$degf)
    if (is.null(d1)) {
      warning("the ", what, " test is not computable: the within-imputation covariance block is ",
              "singular, so the design cannot identify this model. Estimates are still shown, but ",
              "their confidence intervals are not trustworthy.", call. = FALSE)
      return(list(chisq = NA_real_, F = NA_real_, df1 = df1, df2 = degf,
                  p_chisq = NA_real_, p_F = NA_real_, fmi = NA_real_))
    }
    return(list(
      chisq = d1$F * df1,
      F = d1$F,
      df1 = df1,
      df2 = d1$v,
      p_chisq = stats::pchisq(d1$F * df1, df = df1, lower.tail = FALSE),
      p_F = stats::pf(d1$F, df1 = df1, df2 = d1$v, lower.tail = FALSE),
      fmi = d1$r / (1 + d1$r)
    ))
  }

  Vinv <- tryCatch(solve(V), error = function(e) NULL)
  if (is.null(Vinv)) {
    ## Warn rather than stop. Erroring here made svyrcs() refuse a fit that svyrcs_curve() would
    ## happily produce a band for -- the same data giving opposite answers depending on the entry
    ## point. The curve is still meaningful as a point estimate, so report the test as unavailable
    ## and let the caller see the rest.
    warning("the ", what, " test is not computable: the spline covariance block is singular, so ",
            "the design cannot identify this model. This usually means too many knots for the ",
            "amount of information in the data; try fewer knots. Estimates are still shown, but ",
            "their confidence intervals are not trustworthy.", call. = FALSE)
    return(list(chisq = NA_real_, F = NA_real_, df1 = df1, df2 = degf,
                p_chisq = NA_real_, p_F = NA_real_, fmi = if (is.null(mi)) NULL else NA_real_))
  }
  chisq <- as.numeric(t(beta) %*% Vinv %*% beta)
  Fstat <- chisq / df1

  list(
    chisq = chisq,
    F = Fstat,
    df1 = df1,
    df2 = degf,
    p_chisq = stats::pchisq(chisq, df = df1, lower.tail = FALSE),
    p_F = stats::pf(Fstat, df1 = df1, df2 = degf, lower.tail = FALSE),
    fmi = NULL
  )
}

## Overall association and non-linearity tests for a spline block.
##
## The first spline column is x itself, so dropping it leaves exactly the columns that make the fit
## depart from a straight line: testing those jointly is the test of non-linearity.
## Subset an mi bundle to a coefficient block, keeping it in step with the block passed as `V`.
mi_subset <- function(mi, idx) {
  if (is.null(mi)) return(NULL)
  list(Ubar = mi$Ubar[idx, idx, drop = FALSE], B = mi$B[idx, idx, drop = FALSE],
       m = mi$m, degf = mi$degf)
}

## Project an mi bundle through a selection matrix, so it stays in the basis of L V L'.
mi_project <- function(mi, L) {
  if (is.null(mi)) return(NULL)
  list(Ubar = L %*% mi$Ubar %*% t(L), B = L %*% mi$B %*% t(L), m = mi$m, degf = mi$degf)
}

rcs_tests <- function(beta, V, degf, mi = NULL) {
  overall <- wald_block(beta, V, degf, "overall association", mi)
  nonlinear <- if (length(beta) > 1L) {
    wald_block(beta[-1L], V[-1L, -1L, drop = FALSE], degf, "non-linearity", mi_subset(mi, -1L))
  } else {
    list(chisq = NA_real_, F = NA_real_, df1 = 0L, df2 = degf,
         p_chisq = NA_real_, p_F = NA_real_, fmi = if (is.null(mi)) NULL else 0)
  }
  list(overall = overall, nonlinear = nonlinear)
}

## The same two tests for one group, working from the full coefficient vector through that group's
## selection matrix. b_g = L beta and V_g = L V L' carry the main-interaction covariance.
group_tests <- function(beta, V, L, degf, mi = NULL) {
  b <- as.numeric(L %*% beta)
  Vg <- L %*% V %*% t(L)
  rcs_tests(b, Vg, degf, mi_project(mi, L))
}

## Effect-modification tests.
##
## `interaction` asks whether the association differs between groups at all. `shape` drops the
## linear interaction column for each level, so it asks the narrower and usually more interesting
## question: does the *curvature* differ, as opposed to the whole curve being shifted?
interaction_tests <- function(beta, V, modifier, degf, mi = NULL) {
  int_cols <- unlist(modifier$columns, use.names = FALSE)
  shape_cols <- unlist(lapply(modifier$columns, function(cols) cols[-1L]), use.names = FALSE)

  missing <- setdiff(c(int_cols, shape_cols), names(beta))
  if (length(missing)) {
    stop_svyrcs("interaction coefficients missing from the model: ",
                paste(sQuote(missing), collapse = ", "))
  }

  shape <- if (length(shape_cols)) {
    wald_block(beta[shape_cols], V[shape_cols, shape_cols, drop = FALSE], degf,
               "shape interaction", mi_subset(mi, shape_cols))
  } else {
    list(chisq = NA_real_, F = NA_real_, df1 = 0L, df2 = degf,
         p_chisq = NA_real_, p_F = NA_real_, fmi = if (is.null(mi)) NULL else 0)
  }

  list(
    interaction = wald_block(beta[int_cols], V[int_cols, int_cols, drop = FALSE], degf,
                             "interaction", mi_subset(mi, int_cols)),
    shape = shape
  )
}
