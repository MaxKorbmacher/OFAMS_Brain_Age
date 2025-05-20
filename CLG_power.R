# Logistic regression power/sensitivity analysis across models
# Max Korbmacher, 10.12.2024 (with error handling improvements)

# Define data path
datapath <- "/Users/max/Documents/Local/MS/results/"

# Read data
data <- read.csv(paste(datapath, "interrim_data.csv", sep = ""))

# Load required libraries
library(dplyr)
library(datawizard)
library(rlist)
library(WebPower)

# Ensure unique data entries
data <- data_unique(data = data, select = eid)

# 1. POWER ####
# Predictors
predictors <- c("geno", "relapses_12mnths_before_baseline", 
                "NfL..pg.ml.", "edss", "PASAT", "smoking_OFAMS", "BL_BMI", 
                "Omega3_suppl", "BAG_c", "baselineC", "baselineV",
                "Vit_A_0","Vit_D_0", "Vit_E_0",
                "PF", "BP", "VT","MH")

# Generate model formulas up to the 5th order
max_order <- 8
all_combinations <- lapply(1:max_order, function(i) combn(predictors, i, simplify = FALSE))
all_formulas <- unlist(lapply(all_combinations, function(comb_list) {
  lapply(comb_list, function(comb) as.formula(paste("CLG ~ age +Current_DMT+", paste(comb, collapse = " + "))))
}))
# Add optional terms (covariate sex)
all_formulas <- c(all_formulas, lapply(all_formulas, function(f) update(f, . ~ . + sex)))

# Define probability function with error handling
prob <- function(form) {
  tryCatch({
    model <- glm(form, data = data, family = "binomial")
    p0 <- 1 / (1 + exp(-1 * as.numeric(model$coefficients[1])))
    p1 <- 1 / (1 + exp(-1 * sum(model$coefficients)))
    return(c(p0, p1))
  }, error = function(e) {
    message("Error in formula: ", form, " - ", e$message)
    return(c(NA, NA))  # Return NA on error
  })
}

# Calculate probabilities for all formulas
probdf <- lapply(seq_along(all_formulas), function(i) {
  prob(all_formulas[[i]])
})
probdf <- data.frame(do.call(rbind, probdf))
names(probdf) <- c("X1", "X2")  # Rename columns for clarity

# Calculate power with error handling
pwr <- numeric(nrow(probdf))
for (i in 1:nrow(probdf)) {
  tryCatch({
    if (!is.na(probdf$X1[i]) && !is.na(probdf$X2[i])) {
      pwr[i] <- wp.logistic(
        n = 88, 
        p0 = probdf$X1[i], 
        p1 = probdf$X2[i], 
        alpha = 0.05, 
        power = NULL, 
        alternative = "greater", 
        family = "normal"
      )$power
    } else {
      pwr[i] <- NA  # Set NA if probabilities are not valid
    }
  }, error = function(e) {
    message("Error in power calculation for row ", i, " - ", e$message)
    pwr[i] <- NA
  })
}
# add formula description to data
pwr = data.frame(pwr)
pwr$Formula = c(as.character(unlist(all_formulas)))
#
# Write results
write.csv(x = pwr, file = "CLG_power.csv",row.names = F)
print("Power analysis completed.")