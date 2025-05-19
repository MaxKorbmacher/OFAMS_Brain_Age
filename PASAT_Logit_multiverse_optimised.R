# Multiverse for Logistic Regression Models
# ChatGPT aided optimisation, Max Korbmacher 06 December 2024
#
library(dplyr)
library(pscl)
library(parallel)
library(bbmle)
library(pROC)
#
# Load data
data <- read.csv("/Users/max/Documents/Local/MS/results/interrim_data.csv")
#
#
# Predictors
predictors <- c("geno", "relapses_12mnths_before_baseline", 
                "NfL..pg.ml.", "edss", "PASAT", "smoking_OFAMS", "BL_BMI", 
                "Omega3_suppl", "BAG_c", "baselineC", "baselineV",
                "Vit_A_0","Vit_D_0", "Vit_E_0",
                "PF", "BP", "GH", "VT", "SF","MH")
# 
# "Treatment_OFAMS", "RE","RF" >> variability = 0
# Generate model formulas
max_order <- 8
all_combinations <- lapply(1:max_order, function(i) combn(predictors, i, simplify = FALSE))
all_formulas <- unlist(lapply(all_combinations, function(comb_list) {
  lapply(comb_list, function(comb) as.formula(paste("CLG ~ age + Current_DMT +", paste(comb, collapse = " + "))))
}))
# Add optional terms (here, this is the covariate sex)
all_formulas <- c(all_formulas, lapply(all_formulas, function(f) update(f, . ~ . + sex)))
#

# Run models in parallel with error handling
num_cores <- detectCores() - 4

results <- mclapply(all_formulas, function(formula) {
  tryCatch({
    # Fit the model
    model <- glm(formula, family = "binomial", data = data)
    predicted = predict(model,data, type="response")
    
    # Extract coefficients
    coefs <- data.frame(
      Names = names(model$coefficients),
      Beta = model$coefficients,
      SE = summary(model)$coefficients[, 2],
      P = summary(model)$coefficients[, 4],
      Formula = call.to.char(formula),
      AUC = auc(model$data$FLG, predicted)[1]
    )
    
    # Return results
    list(
      formula = as.character(formula),
      pseudoR2 = pscl::pR2(model)["McFadden"],
      AUC = auc(model$data$FLG, predicted)[1],
      coefs = coefs
    )
  }, error = function(e) {
    # Return NULL or some placeholder on error
    message("Error in model: ", formula, " - ", e$message)
    NULL
  })
}, mc.cores = num_cores)

# Filter out NULL results
results <- Filter(Negate(is.null), results)

# Combine results
if (length(results) > 0) {
  model_info <- do.call(rbind, lapply(results, function(res) {
    data.frame(Formula = res$formula, McFaddenR2 = res$pseudoR2, AUC = res$AUC)
  }))
  coef_table <- do.call(rbind, lapply(results, function(res) res$coefs))
  
  # Write results
  write.csv(model_info, "/Users/max/Documents/Local/MS/results/PASAT_model_info_optimised.csv", row.names = FALSE)
  write.csv(coef_table, "/Users/max/Documents/Local/MS/results/PASAT_stratification_optimised.csv", row.names = FALSE)
  
  print("Done.")
} else {
  stop("No models converged or produced results.")
}

# OLD SCRIPT NOT ACCOUNTING FOR ERRORS THROWN BY SINGULARITY
# ---------------------- 
#
# # Run models in parallel
# num_cores <- detectCores() - 4
# results <- mclapply(all_formulas, function(formula) {
#   model <- glm(formula, family = "binomial", data = data)
#   coefs <- data.frame(
#     Names = names(model$coefficients),
#     Beta = model$coefficients,
#     SE = summary(model)$coefficients[, 2],
#     P = summary(model)$coefficients[, 4]
#   )
#   list(
#     formula = as.character(formula),
#     pseudoR2 = pscl::pR2(model)["McFadden"],
#     coefs = coefs
#   )
# }, mc.cores = num_cores)
# 
# # Combine results
# model_info <- do.call(rbind, lapply(results, function(res) {
#   data.frame(Formula = res$formula, McFaddenR2 = res$pseudoR2)
# }))
# coef_table <- do.call(rbind, lapply(results, function(res) res$coefs))
# 
# # Write results
# write.csv(model_info, "/Users/max/Documents/Local/MS/results/EDSS_model_info_optimised.csv", row.names = FALSE)
# write.csv(coef_table, "/Users/max/Documents/Local/MS/results/EDSS_stratification_optimised.csv", row.names = FALSE)
# 
# print("Done.")
