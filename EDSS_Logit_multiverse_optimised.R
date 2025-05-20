# Multiverse for Logistic Regression Models (FLG)
# ChatGPT-aided optimisation – Max Korbmacher, 06 December 2024

library(dplyr)
library(pscl)
library(parallel)
library(pROC)
library(DescTools)

# Load data
data <- read.csv("/Users/max/Documents/Local/MS/results/interrim_data.csv")

# Define predictors
predictors <- c(
  "geno", "relapses_12mnths_before_baseline", 
  "NfL..pg.ml.", "edss", "PASAT", "smoking_OFAMS", "BL_BMI", 
  "Omega3_suppl", "BAG_c", "baselineC", "baselineV",
  "Vit_A_0", "Vit_D_0", "Vit_E_0",
  "PF", "BP", "VT", "MH"
)

# Set maximum number of predictors per model (to avoid overfitting or crashes)
max_order <- 8

# Generate all combinations of predictors up to max_order
all_combinations <- lapply(1:max_order, function(i) combn(predictors, i, simplify = FALSE))

# Generate all base formulas
all_formulas <- unlist(lapply(all_combinations, function(comb_list) {
  lapply(comb_list, function(comb) {
    as.formula(paste("FLG ~ age + Current_DMT +", paste(comb, collapse = " + ")))
  })
}), recursive = FALSE)

# Add optional covariate (sex)
all_formulas <- c(all_formulas, lapply(all_formulas, function(f) update(f, . ~ . + sex)))

# Parallel setup
num_cores <- detectCores() - 4  # leave some cores free

# Run models in parallel with error handling
results <- mclapply(all_formulas, function(formula) {
  tryCatch({
    model <- glm(formula, family = "binomial", data = data)
    predicted <- predict(model, data, type = "response")
    actuals <- data$FLG
    
    # Compute evaluation metrics
    auc_value <- as.numeric(pROC::auc(actuals, predicted))
    brier_value <- DescTools::BrierScore(model)
    
    # Extract coefficients
    coefs <- data.frame(
      Names = names(model$coefficients),
      Beta = model$coefficients,
      SE = summary(model)$coefficients[, 2],
      P = summary(model)$coefficients[, 4],
      Formula = paste(deparse(formula), collapse = ""),
      AUC = auc_value,
      Brier = brier_value
    )
    
    list(
      formula = paste(deparse(formula), collapse = ""),
      pseudoR2 = pscl::pR2(model)["McFadden"],
      AUC = auc_value,
      Brier = brier_value,
      coefs = coefs
    )
  }, error = function(e) {
    message("Error in model: ", paste(deparse(formula), collapse = ""), " - ", e$message)
    NULL
  })
}, mc.cores = num_cores)

# Filter successful models
results <- Filter(Negate(is.null), results)

# Save results if models ran successfully
if (length(results) > 0) {
  
  # Combine model-level summary
  model_info <- do.call(rbind, lapply(results, function(res) {
    data.frame(
      Formula = res$formula,
      McFaddenR2 = res$pseudoR2,
      AUC = res$AUC,
      Brier = res$Brier
    )
  }))
  
  # Combine coefficient-level summary
  coef_table <- do.call(rbind, lapply(results, function(res) res$coefs))
  
  # Write to disk
  write.csv(model_info, "/Users/max/Documents/Local/MS/results/EDSS_model_info_optimised.csv", row.names = FALSE)
  write.csv(coef_table, "/Users/max/Documents/Local/MS/results/EDSS_stratification_optimised.csv", row.names = FALSE)
  
  print("Done.")
} else {
  stop("No models converged or produced results.")
}
