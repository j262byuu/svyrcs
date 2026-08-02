## Locating the effect modifier interacting with the spline term.
##
## Everything here works from a fitted model's `terms` and coefficient names, exactly as
## find_rcs_term() does, so it applies equally to a model svyrcs() built and one the user fitted.

## Column names contributed by a factor's non-reference levels: `sexFemale`, `raceOther`, ...
## Logical variables behave like a factor with levels FALSE/TRUE, so the same construction gives
## `smokerTRUE`.
modifier_levels <- function(x, var) {
  if (is.logical(x)) return(c("FALSE", "TRUE"))
  if (is.character(x)) return(sort(unique(x)))
  if (is.factor(x)) return(levels(x))
  stop_svyrcs("the effect modifier '", var, "' is ", class(x)[1L], ". svyrcs models effect ",
              "modification by groups, so '", var, "' must be a factor, character or logical ",
              "variable. Continuous effect modification is not supported.")
}

#' Find the effect modifier interacting with the spline term
#'
#' Internal helper, documented because its rules determine which formulas `svyrcs()` accepts.
#' Returns `NULL` when the spline term appears on its own, which is what keeps the ungrouped code
#' path identical to a fit without any modifier.
#'
#' @param fit A fitted model containing an [rcs()] term.
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
                "Write rcs(", term$var, ", k) * ", mod, " rather than rcs(", term$var, ", k):", mod,
                ", so that the reference level's curve is estimated.")
  }

  mf <- stats::model.frame(fit)
  modvals <- mf[[mod]]
  if (is.null(modvals)) {
    stop_svyrcs("could not recover the effect modifier '", mod, "' from the fitted model")
  }
  levs <- modifier_levels(modvals, mod)
  if (length(levs) < 2L) {
    stop_svyrcs("the effect modifier '", mod, "' has only one level in the fitted data")
  }

  ## Interaction column names, constructed rather than matched by pattern: the spline term label
  ## contains regex metacharacters, so a pattern built from it is unreliable.
  columns <- list()
  for (lvl in levs[-1L]) {
    tag <- paste0(mod, lvl)
    forward <- paste0(spline_cols, ":", tag)
    reverse <- paste0(tag, ":", spline_cols)
    cols <- if (all(forward %in% coef_names)) {
      forward
    } else if (all(reverse %in% coef_names)) {
      ## `m * rcs(x, k)` puts the modifier first in the coefficient name
      reverse
    } else if (any(forward %in% coef_names) || any(reverse %in% coef_names)) {
      ## Some but not all of the interaction columns were estimated: the fit is rank deficient and
      ## the rest were dropped as aliased. That happens when the spline is collinear with the
      ## modifier within a group -- classically when the modifier is derived from the exposure
      ## itself, or when there are more knots than the smallest group can support.
      present <- intersect(c(forward, reverse), coef_names)
      stop_svyrcs("the model is rank deficient for level '", lvl, "' of '", mod, "': only ",
                  length(present), " of ", length(spline_cols), " interaction coefficients could ",
                  "be estimated, the rest were dropped as aliased. The spline is collinear with '",
                  mod, "' within that group. Use fewer knots, or a modifier that is not derived ",
                  "from '", term$var, "'.")
    } else {
      stop_svyrcs("could not find the interaction coefficients for level '", lvl, "' of '", mod,
                  "'. Expected ", paste(sQuote(forward), collapse = ", "), ".")
    }
    if (anyNA(stats::coef(fit)[cols])) {
      stop_svyrcs("the interaction coefficients for level '", lvl, "' of '", mod, "' include ",
                  "missing values, meaning the fit is rank deficient. Use fewer knots, or a ",
                  "modifier that is not collinear with '", term$var, "' within a group.")
    }
    columns[[lvl]] <- cols
  }

  list(var = mod, levels = levs, ref_level = levs[1L], columns = columns,
       spline_cols = spline_cols)
}

## Map the full coefficient vector to one group's effective spline coefficients.
##
## Row i picks the i-th spline main effect and, for a non-reference level, adds the matching
## interaction coefficient. Building this as a matrix rather than summing two coefficient vectors is
## what makes the variance pick up Cov(main, interaction); adding separate variances would
## understate the standard error for every non-reference group.
selection_matrix <- function(coef_names, spline_cols, int_cols = NULL) {
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
    L[cbind(seq_along(int_cols), inter)] <- 1
  }
  L
}

## The selection matrix for a given group, or the plain spline selection when ungrouped.
group_selection <- function(coef_names, term, modifier, group = NULL) {
  spline_cols <- spline_colnames(term$label, term$nk)
  if (is.null(modifier) || is.null(group) || identical(group, modifier$ref_level)) {
    return(selection_matrix(coef_names, spline_cols, NULL))
  }
  selection_matrix(coef_names, spline_cols, modifier$columns[[group]])
}
