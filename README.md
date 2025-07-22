# Survey-Weighted Restricted Cubic Splines in R (NHANES Tutorial)

## Introduction
Restricted Cubic Splines (RCS) are a great way to model non-linear relationships, especially when simple linear models just don’t cut it. But if you’re working with complex survey data like NHANES, you’ll quickly run into trouble , the usual spline methods don’t know what to do with things like:

1. Sampling weights
2. Clustering
3. Stratification

You need the weights for proper estimation, and the other two are essential for getting valid confidence intervals. 

This little guide walks through my personal workflow for running survey-weighted RCS models with proper variance estimation (yes, Taylor series style). It’s not perfect since I’m not a PhD or anything. So if you spot something off, feel free to open an issue or PR. I’m always learning too.

## Prerequisites
1. R version: 4.0 or higher (even better with Intel MKL)
2. Required packages:
```
r
install.packages(c("survey", "rms", "survival", "tidyverse"))
```

## Why Use RCS?

Real‑world relationships rarely follow a straight line. For example, the link between sedentary behavior and health outcomes is anything but linear. See this excellent [JAMA paper](https://jamanetwork.com/journals/jama/fullarticle/2809418). This paper is simple, solid, and straightforward. I really like it.

### Here’s what the usual approaches do:

1. Categorize continuous variables → Loss of information and power
2. Assume linearity → Model misspecification
3. Use polynomials → Unstable at extremes

### Now here’s what RCS does:

1. Flexible modeling of non-linear relationships
2. Smooth curves with good behavior at extremes
3. Interpretable results (especially with HR/OR)
4. Data-driven knot placement (~~easy to manipulate~~)

## Challenge

Combining RCS with survey design presents unique challenges:

1. No existing integrated solution: The `rms` package doesn't handle survey weights; the `survey` package doesn't have built-in rcs function
2. Reference value complexity: Need to maintain consistent reference values across weighted and unweighted analyses
3. Variance estimation: Must properly account for correlation between predictions
4. Computational burden: Need efficient matrix operations

## My Solution

```
# Survey-weighted RCS with all 3 reference methods
# Including corrected p-value calculations

library(survey)
library(rms)
library(survival)
library(tidyverse)
library(patchwork)

## 0. Setup datadist ---------------------------------------------------------
dataDist <- datadist(select(
  data,
  exposure_var, outcome_time, outcome_event,
  covariate1, covariate2, covariate3  # include all model variables
))
options(datadist = "dataDist")

## 1. Define knots and create spline terms -----------------------------------
knots <- quantile(data$exposure_var, probs = c(.10, .50, .90), na.rm = TRUE)
data$exposure_varRcs <- rcs(data$exposure_var, knots)

## 2. Create survey design ---------------------------------------------------
surveyDesign <- svydesign(
  id      = ~PSU,
  strata  = ~strata,
  weights = ~survey_weight,
  nest    = TRUE,
  data    = data
)

df <- survey::degf(surveyDesign)

## 3. Fit survey weighted Cox model ------------------------------------------
# Example: Cox proportional hazards model
weightedModel <- svycoxph(
  Surv(outcome_time, outcome_event) ~ 
    exposure_rcs + covariate1 + covariate2 + covariate3 + strata(stratification_var),
  design = surveyDesign
)

## 4. Function to calculate p-values correctly -------------------------------
calculate_rcs_pvalues <- function(model, var_name = "exposure_var") {
  
  # Get coefficient names
  coef_names <- names(coef(model))
  
  # Find RCS terms for the variable
  rcs_pattern <- paste0(var_name, "Rcs")
  rcs_terms <- grep(rcs_pattern, coef_names, value = TRUE)
  
  if(length(rcs_terms) == 0) {
    stop("No RCS terms found for ", var_name)
  }
  
  # Get coefficients and variance-covariance matrix
  all_coef <- coef(model)
  vcov_mat <- vcov(model)
  
  # 1. Overall effect (all RCS terms together)
  overall_idx <- which(coef_names %in% rcs_terms)
  overall_coef <- all_coef[overall_idx]
  overall_vcov <- vcov_mat[overall_idx, overall_idx, drop = FALSE]
  
  # Wald test for overall effect
  overall_wald <- as.numeric(t(overall_coef) %*% solve(overall_vcov) %*% overall_coef)
  overall_df <- length(overall_idx)
  overall_pval <- pchisq(overall_wald, df = overall_df, lower.tail = FALSE)
  
  # 2. Non-linear effect (excluding first RCS term)
  if(length(rcs_terms) > 1) {
    nonlin_idx <- overall_idx[-1]  # Exclude first RCS term
    nonlin_coef <- all_coef[nonlin_idx]
    nonlin_vcov <- vcov_mat[nonlin_idx, nonlin_idx, drop = FALSE]
    
    # Wald test for non-linearity
    nonlin_wald <- as.numeric(t(nonlin_coef) %*% solve(nonlin_vcov) %*% nonlin_coef)
    nonlin_df <- length(nonlin_idx)
    nonlin_pval <- pchisq(nonlin_wald, df = nonlin_df, lower.tail = FALSE)
  } else {
    nonlin_wald <- NA
    nonlin_df <- 0
    nonlin_pval <- NA
  }
  
  list(
    overall = list(chisq = overall_wald, df = overall_df, pvalue = overall_pval),
    nonlinear = list(chisq = nonlin_wald, df = nonlin_df, pvalue = nonlin_pval)
  )
}

# Calculate p-values
pvals <- calculate_rcs_pvalues(weightedModel, "exposure_var")

## 5. Create prediction grid -------------------------------------------------
grid_x <- seq(quantile(data$exposure_var, .01, na.rm = TRUE),
              quantile(data$exposure_var, .99, na.rm = TRUE), 
              length.out = 500)

newData <- data.frame(exposure_var = grid_x)
newData$exposure_varRcs <- rcs(newData$exposure_var, knots)

# Set other covariates to reference values
for (v in setdiff(names(data), names(newData))) {
  if (is.factor(data[[v]])) {
    newData[[v]] <- factor(names(sort(table(data[[v]]), decreasing=TRUE)[1]),
                           levels(data[[v]]))
  } else {
    newData[[v]] <- median(data[[v]], na.rm = TRUE)
  }
}

## 6. Function to calculate HR with specific reference -----------------------
calculate_hr <- function(ref_value, ref_name) {
  refData <- newData[1, , drop = FALSE]
  refData$exposure_var <- ref_value
  refData$exposure_varRcs <- rcs(ref_value, knots)
  
  # Design matrices
  tmp <- rbind(newData, refData)
  mm <- model.matrix(weightedModel, data = tmp)
  X <- mm[1:nrow(newData), , drop = FALSE]
  X0 <- mm[nrow(tmp), , drop = FALSE]
  
  # Calculate HR with proper variance
  V <- vcov(weightedModel)
  Xd <- sweep(X, 2, X0)
  var_logHR <- rowSums((Xd %*% V) * Xd)
  se_logHR <- sqrt(var_logHR)
  
  lp_grid <- predict(weightedModel, newdata = newData, type = "lp")
  lp_ref <- as.numeric(predict(weightedModel, newdata = refData, type = "lp"))
  logHR <- lp_grid - lp_ref
  
  crit <- qt(0.975, df)
  HR <- exp(logHR)
  lower_CI <- exp(logHR - crit * se_logHR)
  upper_CI <- exp(logHR + crit * se_logHR)
  
  data.frame(
    exposure_var = newData$exposure_var,
    HR = HR,
    lower = lower_CI,
    upper = upper_CI,
    method = ref_name,
    ref_value = ref_value
  )
}

## 7. Find reference values for each method ----------------------------------

# Method 1: Probability-based (median)
ref_prob <- quantile(data$exposure_var, probs = 0.5, na.rm = TRUE)

# Method 2: Minimum HR as reference 
temp_hr <- calculate_hr(ref_prob, "temp")
ref_minimal <- temp_hr$exposure_var[which.min(temp_hr$HR)]

# Method 3: Maximal HR as reference 
ref_maximal <- temp_hr$exposure_var[which.max(temp_hr$HR)]

## 8. Calculate results for all methods --------------------------------------
results_prob <- calculate_hr(ref_prob, "Probability (Median)")
results_minimal <- calculate_hr(ref_minimal, "Minimal Risk")
results_maximal <- calculate_hr(ref_maximal, "Maximal Risk")

all_results <- rbind(results_prob, results_minimal, results_maximal)

## 9. Create plots for all methods -------------------------------------------
create_plot <- function(data, title_suffix) {
  ref_val <- unique(data$ref_value)
  
  p <- ggplot(data, aes(x = exposure_var, y = HR)) +
    geom_line(linewidth = 1.2, color = "black") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = ref_val, linetype = "dotted", color = "red", alpha = 0.7) +
    geom_point(x = ref_val, y = 1, size = 3, color = "red") +
    scale_x_continuous(breaks = seq(0, 14000, by = 2000)) +
    scale_y_continuous(trans = "log2") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 11)
    ) +
    labs(
      title = paste0(title_suffix, " (Ref: ", round(ref_val, 0), ")"),
      x = "Exposure",
      y = "Hazard Ratio (95% CI)"
    )
  
  return(p)
}

# Create individual plots
p1 <- create_plot(results_prob, "Probability-based (Median)")
p2 <- create_plot(results_minimal, "Minimal Risk")
p3 <- create_plot(results_maximal, "Maximal Risk")

# Combine plots
combined_plot <- (p1 | p2 | p3) +
  plot_annotation(
    title = "Non-linear Association",
    subtitle = paste0("Survey-weighted Cox model with RCS (3 knots); Taylor series 95% CI (df=", df, ")"),
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 12))
  )

print(combined_plot)

## 10. Model Summary ---------------------------------------------------------
cat("\n================== MODEL SUMMARY ==================\n")
cat("Events:", sum(data$outcome_event, na.rm = TRUE), "\n")
cat("Observations:", nrow(data), "\n")
cat("Degrees of freedom:", df, "\n")
cat("\nKnot locations:", round(knots, 0), "steps/day\n")

cat("\n================== P-VALUES =======================\n")
cat("Overall effect of exposure_var:\n")
cat("  Chi-square:", round(pvals$overall$chisq, 2), "\n")
cat("  d.f.:", pvals$overall$df, "\n")
cat("  P-value:", format.pval(pvals$overall$pvalue, digits = 4), "\n")

cat("\nNon-linear component:\n")
cat("  Chi-square:", round(pvals$nonlinear$chisq, 2), "\n")
cat("  d.f.:", pvals$nonlinear$df, "\n")
cat("  P-value:", format.pval(pvals$nonlinear$pvalue, digits = 4), "\n")

cat("\n================== REFERENCE VALUES ===============\n")
cat("Probability (median):", round(ref_prob, 0), "steps/day\n")
cat("Minimal risk:", round(ref_minimal, 0), "steps/day\n")
cat("Maximal risk:", round(ref_maximal, 0), "steps/day\n")

cat("\n================== LOWEST/HIGHEST POINTS ==========\n")
# Find global minimum and maximum HRs across all methods
global_results <- calculate_hr(ref_prob, "global")
min_idx <- which.min(global_results$HR)
max_idx <- which.max(global_results$HR)

cat("Lowest HR point:", round(global_results$exposure_var[min_idx], 0), 
    "steps/day (HR =", round(global_results$HR[min_idx], 3), ")\n")
cat("Highest HR point:", round(global_results$exposure_var[max_idx], 0), 
    "steps/day (HR =", round(global_results$HR[max_idx], 3), ")\n")

cat("\n================== INTERPRETATION =================\n")
if(pvals$overall$pvalue < 0.05) {
  cat("exposure_var is significantly associated with mortality (p =", 
      format.pval(pvals$overall$pvalue, digits = 4), ")\n")
} else {
  cat("exposure_var is NOT significantly associated with mortality (p =", 
      format.pval(pvals$overall$pvalue, digits = 4), ")\n")
}

if(!is.na(pvals$nonlinear$pvalue) && pvals$nonlinear$pvalue < 0.05) {
  cat("The relationship is significantly non-linear (p =", 
      format.pval(pvals$nonlinear$pvalue, digits = 4), ")\n")
  cat("A flexible model (RCS) is justified over a simple linear model\n")
} else {
  cat("The relationship is NOT significantly non-linear (p =", 
      format.pval(pvals$nonlinear$pvalue, digits = 4), ")\n")
  cat("A linear model might be sufficient\n")
}
```
## Common Troubleshooting
1. Adapting for Logistic Regression

```
r
# For binary outcomes:
weightedModel <- svyglm(
  outcome ~ exposure_rcs + covariate1 + covariate2 + covariate3,
  design = surveyDesign,
  family = quasibinomial()  # For proper variance estimation
)

# Then change HR to OR in the results:
OR <- exp(logOR)
# Everything else stays the same!
```

2. Different Number of Knots

```
r
# For 5 knots (more flexibility):
knots <- quantile(data$exposure_var, 
                  probs = c(0.05, 0.275, 0.50, 0.725, 0.95), 
                  na.rm = TRUE)

# For 4 knots:
knots <- quantile(data$exposure_var, 
                  probs = c(0.05, 0.35, 0.65, 0.95), 
                  na.rm = TRUE)

# Remember: more knots = more flexibility but less stability
```
3. Missing Data

```
# Option 1: Complete case analysis (default)
surveyDesign <- svydesign(..., data = na.omit(data))

# Option 2: Multiple imputation (recommended)
# See the survey package vignette on svymi
```

## Citations
If you used this guide's svyrcs approach, I would appreciate if you cite this page! Thanks!

## References
1. **Restricted Cubic Splines**: Harrell FE Jr. *Regression Modeling Strategies*. 2nd ed. Springer; 2015. [doi:10.1007/978-3-319-19425-7](https://doi.org/10.1007/978-3-319-19425-7)
2. **Survey Package**: Lumley T. *Complex Surveys: A Guide to Analysis Using R*. Wiley; 2010. [doi:10.1002/9780470580066](https://doi.org/10.1002/9780470580066)
3. **Taylor Series Linearization**: Binder DA. On the variances of asymptotically normal estimators from complex surveys. *International Statistical Review*. 1983;51(3):279-292. [doi:10.2307/1402588](https://doi.org/10.2307/1402588)
4. **NHANES Analysis Guidelines**: [CDC NHANES Tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/default.aspx)

## Contributing: 
Found an error? Have a better approach? Please open an issue or submit a PR!
