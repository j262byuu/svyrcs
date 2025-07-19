# My Personal Guide to Using Restricted Cubic Splines in Complex Survey Data (a.k.a. NHANES)

## Introduction
Restricted Cubic Splines (RCS) are a great way to model non-linear relationships, especially when simple linear models just don’t cut it. But if you’re working with complex survey data like NHANES, you’ll quickly run into trouble , the usual spline methods don’t know what to do with things like:

1. Sampling weights
2. Clustering
3. Stratification

You need the weights for proper estimation, and the other two are essential for getting valid confidence intervals. 

This little guide walks through my personal workflow for running survey-weighted RCS models with proper variance estimation (yes, Taylor series style). It’s not perfect since I’m not a PhD or anything. So if you spot something off, feel free to open an issue or PR. I’m always learning too.

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

1. No existing integrated solution: The rms package doesn't handle survey weights; the survey package doesn't have built-in rcs function
2. Reference value complexity: Need to maintain consistent reference values across weighted and unweighted analyses
3. Variance estimation: Must properly account for correlation between predictions
4. Computational burden: Need efficient matrix operations

## My Solution

```
# -------------------------------------------------------------------------
# Survey weighted restricted cubic spline with proper variance estimation
# Using Taylor linearization with contrast matrices
# -------------------------------------------------------------------------

library(survey)
library(rms)
library(survival)
library(tidyverse)

## 1. Setup data distribution for reference values ------------------------
dataDist <- datadist(select(
  data,
  exposure_var, outcome_time, outcome_event,
  covariate1, covariate2, covariate3  # include all model variables
))
options(datadist = "dataDist")

## 2. Define knots for RCS ------------------------------------------------
# Using Harrell's recommended percentiles for 3 knots
knots <- quantile(data$exposure_var, probs = c(0.10, 0.50, 0.90), na.rm = TRUE)

## 3. Create spline terms in the dataset ----------------------------------
data$exposure_rcs <- rcs(data$exposure_var, knots)

## 4. Create survey design object -----------------------------------------
surveyDesign <- svydesign(
  id      = ~PSU,
  strata  = ~strata,
  weights = ~survey_weight,
  nest    = TRUE,
  data    = data
)

# Get design degrees of freedom for proper inference
df <- survey::degf(surveyDesign)

## 5. Fit survey-weighted model -------------------------------------------
# Example: Cox proportional hazards model
weightedModel <- svycoxph(
  Surv(outcome_time, outcome_event) ~ 
    exposure_rcs + covariate1 + covariate2 + covariate3 + strata(stratification_var),
  design = surveyDesign
)

## 6. Create prediction grid ----------------------------------------------
# Grid from 1st to 99th percentiles (avoid sparse extremes)
grid_values <- seq(
  quantile(data$exposure_var, 0.01, na.rm = TRUE),
  quantile(data$exposure_var, 0.99, na.rm = TRUE), 
  length.out = 500
)

# Create prediction dataset
newData <- data.frame(exposure_var = grid_values)
newData$exposure_rcs <- rcs(newData$exposure_var, knots)

# Set other covariates to reference values
for (v in setdiff(all.vars(formula(weightedModel)), names(newData))) {
  if (v %in% names(data)) {
    if (is.factor(data[[v]])) {
      # Use mode for factors
      tbl <- table(data[[v]], useNA = "no")
      newData[[v]] <- factor(
        names(sort(tbl, decreasing = TRUE)[1]),
        levels = levels(data[[v]])
      )
    } else if (is.numeric(data[[v]])) {
      # Use median for continuous
      newData[[v]] <- median(data[[v]], na.rm = TRUE)
    }
  }
}

## 7. Define reference exposure value -------------------------------------
# Get reference from datadist (typically median)
refValue <- dataDist$limits["Adjust to", "exposure_var"]
refData <- newData[1, , drop = FALSE]
refData$exposure_var <- refValue
refData$exposure_rcs <- rcs(refValue, knots)

## 8. Calculate predictions using contrast matrices -----------------------
# This is the KEY step for proper variance estimation

# Combine prediction and reference data
tmp <- rbind(newData, refData)

# Get model matrices
mm <- model.matrix(weightedModel, data = tmp)
X <- mm[1:nrow(newData), , drop = FALSE]    # Prediction points
X0 <- mm[nrow(tmp), , drop = FALSE]         # Reference point

# Get variance-covariance matrix from survey-weighted model
V <- vcov(weightedModel)

# Calculate contrasts: X - X0
Xd <- sweep(X, 2, X0)

# Variance of log hazard ratios: (X-X0)' V (X-X0)
# This properly accounts for correlation!
var_logHR <- rowSums((Xd %*% V) * Xd)
se_logHR <- sqrt(var_logHR)

# Get linear predictors
lp_grid <- predict(weightedModel, newdata = newData, type = "lp")
lp_ref <- as.numeric(predict(weightedModel, newdata = refData, type = "lp"))
logHR <- lp_grid - lp_ref

## 9. Calculate confidence intervals --------------------------------------
# Use t-distribution with survey design degrees of freedom
crit <- qt(0.975, df)

# Transform to hazard ratio scale
HR <- exp(logHR)
lower_CI <- exp(logHR - crit * se_logHR)
upper_CI <- exp(logHR + crit * se_logHR)

# Compile results
results <- data.frame(
  exposure = newData$exposure_var,
  HR = HR,
  lower = lower_CI,
  upper = upper_CI,
  se_logHR = se_logHR
)

## 10. Visualize results --------------------------------------------------
ggplot(results, aes(x = exposure, y = HR)) +
  geom_line(size = 1.2, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  # Use log scale for y-axis (recommended for ratios)
  scale_y_continuous(
    trans = "log2",
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = c(0.25, 0.5, 1, 2, 4, 8)
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  labs(
    x = "Exposure",
    y = "Hazard Ratio (95% CI)",
    title = "Non-linear Association",
    subtitle = sprintf("Survey-weighted Cox model with RCS (3 knots); Taylor series 95%% CI (df=%d)", df)
  )

## 11. Report key statistics ----------------------------------------------
# Hazard ratios at selected exposure values
key_exposures <- c(10, 25, 50, 75, 90)  # percentiles or meaningful values
key_indices <- sapply(key_exposures, function(x) which.min(abs(results$exposure - x)))

summary_table <- results[key_indices, ]
summary_table$CI_width <- summary_table$upper - summary_table$lower
print(round(summary_table, 3))

# Model information
cat("\nModel Summary:\n")
cat("Events:", sum(data$outcome_event, na.rm = TRUE), "\n")
cat("Observations:", nrow(data), "\n")
cat("Reference value:", round(refValue, 2), "\n")
cat("Knot locations:", round(knots, 2), "\n")
cat("Degrees of freedom:", df, "\n")
```
