# Load required packages
library(tidyverse)
library(car)
library(fixest)
library(terra)

# Load processed data (new standardized path)
data_all <- qread("data/derived/tidied.qs", nthreads = qn) %>% drop_na()

# Filter for Rice crops (both LR and SR&ER)
adata <- data_all %>% filter(grepl("Rice", crop_parent))

# Define base regression formula consistent with feols
fml_base <- as.formula(paste(
  "GOSIF_sum ~",
  "Caod + I(Caod^2) +",
  "Faod + I(Faod^2) +",
  "cloud + I(cloud^2) +",
  "AOT40 + Faod:AOT40 + I(Faod^2):AOT40 +",
  paste0("bin", 1:42, collapse = " + "), "+",
  paste0("root_", 1:9, collapse = " + "), "+",
  paste0("d", 1:13, collapse = " + ")
))

# Build model matrix
X <- model.matrix(fml_base, data = adata)
X_noint <- X[, !colnames(X) %in% "(Intercept)"]

# Identify alias (perfectly collinear) variables
mod_tmp <- lm(GOSIF_sum ~ . - 1, data = cbind(GOSIF_sum = adata$GOSIF_sum, X_noint))
alias_vars <- names(alias(mod_tmp)$Complete)

# Remove completely collinear variables
X_clean <- X_noint[, !colnames(X_noint) %in% alias_vars]

# Fit clean linear model and compute VIF
mod_lm_clean <- lm(GOSIF_sum ~ . - 1, data = cbind(GOSIF_sum = adata$GOSIF_sum, X_clean))
vif_vals <- vif(mod_lm_clean)

cat("❗ Variables removed due to complete collinearity:\n")
print(alias_vars)
cat("\nVIF values for remaining variables:\n")
print(vif_vals)

# Focus VIF on Faod, AOT40, and interaction terms
X_focus <- X_clean %>% select(matches("^Faod$|^AOT40$|Faod:AOT40"))

# Compute focused VIF
mod_focus <- lm(GOSIF_sum ~ . - 1, data = cbind(GOSIF_sum = adata$GOSIF_sum, X_focus))
vif_focus <- vif(mod_focus)

cat("\nVIF values for Faod, AOT40 and interactions:\n")
print(vif_focus)

# Save results to new standardized directory
write.csv(as.data.frame(vif_vals), "data/derived/vif_results_rice.csv", row.names = TRUE)
write.csv(as.data.frame(vif_focus), "data/derived/vif_focus_rice.csv", row.names = TRUE)