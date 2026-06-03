# Load packages and helper functions
source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

library(dplyr)
library(fixest)
library(purrr)
library(qs)
library(tidyverse)

# Define paths
data_dir   <- "data/derived/crop/maize"
output_dir <- "data/derived/ablation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load final panel data
fnl_data <- qread(file.path(data_dir, "tidied.qs"))

# Define modules including CO2
module_defs <- list(
  Caod  = c("Caod", "Caod^2"),
  Faod  = c("Faod", "Faod^2"),
  AOT40 = "AOT40",
  cloud = c("cloud", "cloud^2"),
  bins  = str_c("bin", 1:42),
  roots = str_c("root_", 1:9),
  ds    = str_c("d", 1:13),
  co2   = "co2",  # Added CO2 as a module
  interact = c("Faod * AOT40", "(Faod^2) * AOT40")
)

# Generate all possible module inclusion/exclusion combinations
module_grid <- expand.grid(
  Caod     = c(TRUE, FALSE),
  Faod     = c(TRUE, FALSE),
  AOT40    = c(TRUE, FALSE),
  cloud    = c(TRUE, FALSE),
  bins     = c(TRUE, FALSE),
  roots    = c(TRUE, FALSE),
  ds       = c(TRUE, FALSE),
  co2      = c(TRUE, FALSE),  # Include CO2
  interact = c(TRUE, FALSE)
)

# Remove duplicate or invalid combinations if needed
module_grid <- module_grid %>% distinct()

# Function to build formula dynamically based on modules
build_formula <- function(modules, module_defs, response = "GOSIF_sum") {
  included_vars <- unlist(module_defs[names(modules)[modules]])
  formula_str <- paste(response, "~", paste(included_vars, collapse = "+"), "| x_y[year]")
  as.formula(formula_str)
}

# Loop over all module combinations and run regression
results_list <- list()

for (i in seq_len(nrow(module_grid))) {
  
  modules <- module_grid[i, ]
  formula_i <- build_formula(modules, module_defs)
  
  model_i <- feols(
    formula_i,
    data = fnl_data,
    weights = ~fraction,
    nthreads = 0,
    notes = FALSE,
    lean = TRUE
  )
  
  # Extract coefficients
  coef_i <- tidy(model_i) %>%
    mutate(model_id = i)
  
  # Extract within R^2
  r2_i <- data.frame(
    model_id = i,
    within_r2 = glance(model_i)$within.r.squared
  )
  
  # Combine
  results_list[[i]] <- list(coef = coef_i, r2 = r2_i, modules = modules)
}

# Combine all results
coef_results <- map_dfr(results_list, "coef")
r2_results <- map_dfr(results_list, "r2")
modules_table <- map_dfr(results_list, "modules") %>%
  mutate(model_id = row_number())

# Merge module info with R²
ablation_results <- modules_table %>%
  left_join(r2_results, by = "model_id") %>%
  arrange(desc(within_r2))

# Save outputs
write.csv(coef_results, file.path(output_dir, "ablation_coefficients.csv"), row.names = FALSE)
write.csv(ablation_results, file.path(output_dir, "ablation_results.csv"), row.names = FALSE)