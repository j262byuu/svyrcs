## Locating the effect modifier interacting with the spline term.
##
## Everything here works from a fitted model's `terms` and coefficient names, exactly as
## find_rcs_term() does, so it applies equally to a model svyrcs() built and one the user fitted.

## The modifier's levels and its contrast matrix.
##
## Working from the contrast matrix rather than the level names is what makes non-default contrasts
## work. Treatment contrasts (the default for an unordered factor) name their columns after the
## non-reference levels, so `paste0(var, colnames(C))` gives `sexFemale`; but an *ordered* factor
## uses polynomial contrasts whose columns are `.L`, `.Q`, `.C`, and a level name will never appear
## in a coefficient name. Level j's effective spline coefficients are then
## `beta_main + sum_k C[j, k] * beta_interaction_k`, which reduces to the 0/1 case for treatment
## contrasts.
modifier_contrasts <- function(x, var, fitted = NULL) {
  if (is.logical(x)) x <- factor(x, levels = c("FALSE", "TRUE"))
  if (is.character(x)) x <- factor(x)
  if (!is.factor(x)) {
    stop_svyrcs("the effect modifier '", var, "' is ", class(x)[1L], ". svyrcs models effect ",
                "modification by groups, so '", var, "' must be a factor, character or logical ",
                "variable. Continuous effect modification is not supported.")
  }
  levs <- levels(x)
  C <- if (!is.null(fitted) && is.matrix(fitted) && nrow(fitted) == length(levs)) {
    fitted
  } else {
    stats::contrasts(x)
  }
  if (is.null(dim(C))) C <- matrix(C, nrow = length(levs))
  rownames(C) <- levs
  if (is.null(colnames(C))) colnames(C) <- as.character(seq_len(ncol(C)))
  list(levels = levs, contrasts = C)
}

#' Find the effect modifier interacting with the spline term
#'
#' Internal helper, documented because its rules determine which formulas `svyrcs()` accepts.
#' Returns `NULL` when the spline term appears on its own, which is what keeps the ungrouped code
#' path identical to a fit without any modifier.
#'
#' @param fit A fitted model containing an [rcspline()] term.
#' @param term The term description from `find_rcs_term()`.
#'
#' @return `NULL`, or a list with `var`, `levels`, `ref_level` and `columns` (a named list of
#'   interaction column names, one entry per non-reference level).
#'
#' @keywords internal
#' @noRd
find_modifier <- function(fit, term) {
  tt <- stats::terms(fit)
  factors <- attr(tt, "factors")
  order <- attr(tt, "order")
  coef_names <- names(stats::coef(fit))

  if (is.null(factors) || !length(factors)) return(NULL)
  row <- match(term$label, rownames(factors))
  if (is.na(row)) return(NULL)

  ## Interaction terms that involve the spline.
  involved <- which(factors[row, ] != 0 & order > 1L)
  if (!length(involved)) return(NULL)
  if (length(involved) > 1L) {
    stop_svyrcs("the spline term for '", term$var, "' is interacted with more than one term (",
                paste(colnames(factors)[involved], collapse = ", "),
                "). svyrcs supports one effect modifier at a time.")
  }

  col <- involved
  partners <- setdiff(rownames(factors)[factors[, col] != 0], term$label)
  if (length(partners) != 1L) {
    stop_svyrcs("the interaction '", colnames(factors)[col], "' involves ", length(partners),
                " variables besides the spline. svyrcs supports a two-way interaction with a ",
                "single effect modifier.")
  }
  mod <- partners

  ## The main effects must be present, otherwise the reference level's curve is misspecified.
  spline_cols <- spline_colnames(term$label, term$nk)
  if (!all(spline_cols %in% coef_names)) {
    stop_svyrcs("the model has an interaction with '", mod, "' but no main effect for the spline. ",
                "Write rcspline(", term$var, ", k) * ", mod, " rather than rcspline(", term$var, ", k):", mod,
                ", so that the reference level's curve is estimated.")
  }

  mf <- stats::model.frame(fit)
  modvals <- mf[[mod]]
  if (is.null(modvals)) {
    stop_svyrcs("could not recover the effect modifier '", mod, "' from the fitted model")
  }
  ## fit$contrasts records the coding the model matrix was actually built with, including one
  ## supplied through contrasts.arg rather than attached to the factor. Fall back to the factor's own
  ## contrasts when the fit does not keep them.
  fitted_contr <- tryCatch(fit$contrasts[[mod]], error = function(e) NULL)
  mc <- modifier_contrasts(modvals, mod, fitted_contr)
  levs <- mc$levels
  C <- mc$contrasts
  if (length(levs) < 2L) {
    stop_svyrcs("the effect modifier '", mod, "' has only one level in the fitted data")
  }

  ## Interaction column names, constructed rather than matched by pattern: the spline term label
  ## contains regex metacharacters, so a pattern built from it is unreliable.
  ## One block of interaction columns per *contrast* column, not per level.
  columns <- list()
  for (cn in colnames(C)) {
    tag <- paste0(mod, cn)
    forward <- paste0(spline_cols, ":", tag)
    reverse <- paste0(tag, ":", spline_cols)
    cols <- if (all(forward %in% coef_names)) {
      forward
    } else if (all(reverse %in% coef_names)) {
      ## `m * rcspline(x, k)` puts the modifier first in the coefficient name
      reverse
    } else if (any(forward %in% coef_names) || any(reverse %in% coef_names)) {
      ## Some but not all of the interaction columns were estimated: the fit is rank deficient and
      ## the rest were dropped as aliased. That happens when the spline is collinear with the
      ## modifier within a group -- classically when the modifier is derived from the exposure
      ## itself, or when there are more knots than the smallest group can support.
      present <- intersect(c(forward, reverse), coef_names)
      stop_svyrcs("the model is rank deficient for '", tag, "': only ", length(present), " of ",
                  length(spline_cols), " interaction coefficients could be estimated, the rest ",
                  "were dropped as aliased. The spline is collinear with '", mod, "' within a ",
                  "group. Use fewer knots, or a modifier that is not derived from '", term$var,
                  "'.")
    } else {
      stop_svyrcs("could not find the interaction coefficients for '", tag, "' of '", mod, "'. ",
                  "Expected ", paste(sQuote(forward), collapse = ", "), ".")
    }
    if (anyNA(stats::coef(fit)[cols])) {
      stop_svyrcs("the interaction coefficients for '", tag, "' of '", mod, "' include missing ",
                  "values, meaning the fit is rank deficient. Use fewer knots, or a modifier that ",
                  "is not collinear with '", term$var, "' within a group.")
    }
    columns[[cn]] <- cols
  }

  list(var = mod, levels = levs, ref_level = levs[1L], columns = columns,
       contrasts = C, spline_cols = spline_cols)
}

## Map the full coefficient vector to one group's effective spline coefficients.
##
## Row i picks the i-th spline main effect and, for a non-reference level, adds the matching
## interaction coefficient. Building this as a matrix rather than summing two coefficient vectors is
## what makes the variance pick up Cov(main, interaction); adding separate variances would
## understate the standard error for every non-reference group.
selection_matrix <- function(coef_names, spline_cols, int_cols = NULL, weight = 1) {
  L <- matrix(0, nrow = length(spline_cols), ncol = length(coef_names),
              dimnames = list(spline_cols, coef_names))
  main <- match(spline_cols, coef_names)
  if (anyNA(main)) {
    stop_svyrcs("spline coefficients missing from the model: ",
                paste(sQuote(spline_cols[is.na(main)]), collapse = ", "))
  }
  L[cbind(seq_along(spline_cols), main)] <- 1
  if (!is.null(int_cols)) {
    inter <- match(int_cols, coef_names)
    if (anyNA(inter)) {
      stop_svyrcs("interaction coefficients missing from the model: ",
                  paste(sQuote(int_cols[is.na(inter)]), collapse = ", "))
    }
    L[cbind(seq_along(int_cols), inter)] <- weight
  }
  L
}

## The selection matrix for a given group, or the plain spline selection when ungrouped.
##
## Level j contributes weight C[j, k] to interaction block k. Under treatment contrasts that is 0 for
## the reference level and a single 1 elsewhere, so this reduces exactly to the previous behaviour;
## under polynomial contrasts (ordered factors) every level draws on every block.
group_selection <- function(coef_names, term, modifier, group = NULL) {
  spline_cols <- spline_colnames(term$label, term$nk)
  L <- selection_matrix(coef_names, spline_cols, NULL)
  if (is.null(modifier) || is.null(group)) return(L)

  w <- modifier$contrasts[group, ]
  for (k in seq_along(modifier$columns)) {
    if (w[[k]] == 0) next
    idx <- match(modifier$columns[[k]], coef_names)
    if (anyNA(idx)) {
      stop_svyrcs("interaction coefficients missing from the model: ",
                  paste(sQuote(modifier$columns[[k]][is.na(idx)]), collapse = ", "))
    }
    L[cbind(seq_along(idx), idx)] <- w[[k]]
  }
  L
}
