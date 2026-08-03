## Multiply imputed fixtures.
##
## The imputations are built deterministically rather than with mice: the suite has to be fast and
## exactly reproducible, and what is being tested is the pooling, not the imputation model. Each
## imputation fills the missing tchol / hdl values by sampling from the observed ones under its own
## fixed seed, which gives genuine between-imputation variance without any modelling.

deterministic_imputations <- function(m = 5, seed = 1) {
  d <- svyrcs::nhanes_bmi
  lapply(seq_len(m), function(i) {
    set.seed(seed + i)
    out <- d
    for (v in c("tchol", "hdl")) {
      miss <- is.na(out[[v]])
      if (any(miss)) out[[v]][miss] <- sample(out[[v]][!miss], sum(miss), replace = TRUE)
    }
    out
  })
}

## A svyimputationList whose model variables genuinely differ between imputations.
mi_design <- local({
  cache <- list()
  function(m = 5, seed = 1) {
    key <- paste(m, seed)
    if (is.null(cache[[key]])) {
      il <- mitools::imputationList(deterministic_imputations(m, seed))
      cache[[key]] <<- survey::svydesign(id = ~psu, strata = ~strata, weights = ~weight,
                                         nest = TRUE, data = il)
    }
    cache[[key]]
  }
})

## A svyimputationList with identical components: nothing is imputed, so the pooled answer must
## reproduce the complete-data fit exactly.
mi_design_degenerate <- local({
  cached <- NULL
  function(m = 3) {
    if (is.null(cached)) {
      il <- mitools::imputationList(rep(list(svyrcs::nhanes_bmi), m))
      cached <<- survey::svydesign(id = ~psu, strata = ~strata, weights = ~weight,
                                   nest = TRUE, data = il)
    }
    cached
  }
})

## Fit the same model in every component design, for tests that pool by hand.
component_fits <- function(designs, formula, family = NULL, cox = FALSE) {
  lapply(designs$designs, function(d) {
    if (cox) survey::svycoxph(formula, design = d)
    else survey::svyglm(formula, design = d, family = family %||% stats::gaussian())
  })
}
