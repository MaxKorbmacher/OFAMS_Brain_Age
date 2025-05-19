# Logistic regression power/sensitivity analysis with formulas in output

# Define data path
datapath <- "/Users/max/Documents/Local/MS/results/"
data <- read.csv(file.path(datapath, "interrim_data.csv"))

# Packages
library(dplyr)
library(datawizard)
library(rlist)
#install.packages("WebPower")
library(WebPower)
library(parallel)

# Ensure unique entries
data <- data_unique(data = data, select = eid)

# Predictors
predictors <- c(
  "geno", "relapses_12mnths_before_baseline", "CH3L.1..mg.ml..mean", 
  "NfL..pg.ml.", "edss", "PASAT", "smoking_OFAMS", "BL_BMI", 
  "Omega3_suppl", "BAG_c", "baselineC", "baselineV", 
  "PF", "RF", "BP", "GH", "VT", "SF", "MH", "Vit_A_0", "Vit_D_0", 
  "Vit_E_0"
) # "RE" >> variability = 0

# Generate model formulas up to the 8th order
max_order <- 8
all_combinations <- unlist(lapply(1:max_order, function(i) combn(predictors, i, simplify = FALSE)), recursive = FALSE)
all_formulas <- lapply(all_combinations, function(comb) {
  as.formula(paste("FLG ~ age + Current_DMT +", paste(comb, collapse = " + ")))
})
# Add optional term (covariate sex)
all_formulas <- c(all_formulas, lapply(all_formulas, function(f) update(f, . ~ . + sex)))

# Function to calculate probabilities
calculate_prob <- function(formula, data) {
  tryCatch({
    model <- glm(formula, data = data, family = "binomial")
    p0 <- 1 / (1 + exp(-model$coefficients[1]))
    p1 <- 1 / (1 + exp(-sum(model$coefficients)))
    return(c(p0 = p0, p1 = p1))
  }, error = function(e) {
    return(c(p0 = NA, p1 = NA)) # Return NA if the model fails
  })
}

# Function to calculate power
calculate_power <- function(probabilities) {
  tryCatch({
    if (any(is.na(probabilities))) return(NA) # Skip invalid probabilities
    wp.logistic(
      n = 85, 
      p0 = probabilities[1], 
      p1 = probabilities[2], 
      alpha = 0.05,
      power = NULL, 
      alternative = "greater", 
      family = "normal"
    )$power
  }, error = function(e) {
    return(NA) # Return NA if power calculation fails
  })
}

# Parallel processing with formulas included in output
cl <- makeCluster(detectCores() - 1) # Use all but one core
clusterExport(cl, c("data", "calculate_prob", "calculate_power", "glm", "wp.logistic"))
clusterEvalQ(cl, library(WebPower))

results <- parLapply(cl, seq_along(all_formulas), function(i) {
  formula <- all_formulas[[i]]
  probs <- calculate_prob(formula, data)
  power <- calculate_power(probs)
  return(list(
    Formula = as.character(formula), # Convert formula to string
    Power = power
  ))
})

stopCluster(cl)

# Convert results to a data frame
results_df <- do.call(rbind, lapply(results, as.data.frame))
results_df <- results_df[!is.na(results_df$Power), ] # Filter valid results

# Save results to CSV
write.csv(x = results_df, file = "power.csv", row.names = FALSE)
