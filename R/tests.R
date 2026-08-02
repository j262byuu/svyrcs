## Joint Wald test on a coefficient block.
##
## Reported both as a chi-square and as a design-based F. With a survey design the chi-square
## reference distribution ignores the fact that the variance is estimated from a modest number of
## PSUs, so it is anti-conservative; F on (df1, degf) is the survey convention and reproduces
## survey::regTermTest() exactly.
wald_block <- function(beta, V, degf, what) {
  df1 <- length(beta)
  if (df1 < 1L) {
    return(list(chisq = NA_real_, F = NA_real_, df1 = 0L, df2 = degf,
                p_chisq = NA_real_, p_F = NA_real_))
  }
  Vinv <- tryCatch(solve(V), error = function(e) NULL)
  if (is.null(Vinv)) {
    stop_svyrcs("the ", what, " test is not computable: the spline covariance block is singular. ",
                "This usually means too many knots for the amount of information in the data; ",
                "try fewer knots.")
  }
  chisq <- as.numeric(t(beta) %*% Vinv %*% beta)
  Fstat <- chisq / df1
  list(
    chisq = chisq,
    F = Fstat,
    df1 = df1,
    df2 = degf,
    p_chisq = stats::pchisq(chisq, df = df1, lower.tail = FALSE),
    p_F = stats::pf(Fstat, df1 = df1, df2 = degf, lower.tail = FALSE)
  )
}

## Overall association and non-linearity tests for a spline block.
##
## The first spline column is x itself, so dropping it leaves exactly the columns that make the fit
## depart from a straight line: testing those jointly is the test of non-linearity.
rcs_tests <- function(beta, V, degf) {
  overall <- wald_block(beta, V, degf, "overall association")
  nonlinear <- if (length(beta) > 1L) {
    wald_block(beta[-1L], V[-1L, -1L, drop = FALSE], degf, "non-linearity")
  } else {
    list(chisq = NA_real_, F = NA_real_, df1 = 0L, df2 = degf,
         p_chisq = NA_real_, p_F = NA_real_)
  }
  list(overall = overall, nonlinear = nonlinear)
}

## The same two tests for one group, working from the full coefficient vector through that group's
## selection matrix. b_g = L beta and V_g = L V L' carry the main-interaction covariance.
group_tests <- function(beta, V, L, degf) {
  b <- as.numeric(L %*% beta)
  Vg <- L %*% V %*% t(L)
  rcs_tests(b, Vg, degf)
}

## Effect-modification tests.
##
## `interaction` asks whether the association differs between groups at all. `shape` drops the
## linear interaction column for each level, so it asks the narrower and usually more interesting
## question: does the *curvature* differ, as opposed to the whole curve being shifted?
interaction_tests <- function(beta, V, modifier, degf) {
  int_cols <- unlist(modifier$columns, use.names = FALSE)
  shape_cols <- unlist(lapply(modifier$columns, function(cols) cols[-1L]), use.names = FALSE)

  missing <- setdiff(c(int_cols, shape_cols), names(beta))
  if (length(missing)) {
    stop_svyrcs("interaction coefficients missing from the model: ",
                paste(sQuote(missing), collapse = ", "))
  }

  shape <- if (length(shape_cols)) {
    wald_block(beta[shape_cols], V[shape_cols, shape_cols, drop = FALSE], degf,
               "shape interaction")
  } else {
    list(chisq = NA_real_, F = NA_real_, df1 = 0L, df2 = degf,
         p_chisq = NA_real_, p_F = NA_real_)
  }

  list(
    interaction = wald_block(beta[int_cols], V[int_cols, int_cols, drop = FALSE], degf,
                             "interaction"),
    shape = shape
  )
}
