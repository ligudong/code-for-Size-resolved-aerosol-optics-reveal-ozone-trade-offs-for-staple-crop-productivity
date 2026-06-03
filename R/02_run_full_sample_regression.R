# Load packages and helper functions
source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

library(dplyr)
library(fixest)
library(broom)
library(purrr)
library(qs)

# Directories
data_dir   <- "data/derived/crop/maize"
output_dir <- "data/derived/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

options(scipen = 999)  # Disable scientific notation

# Load final panel data
fnl_data <- qread(file.path(data_dir, "tidied.qs"))

# Function to run fixed-effect regression per crop
run_full_sample_regression <- function(data, formula_object) {
  
  data %>%
    nest(fdata = -crop_parent) %>%
    mutate(
      coefs = map(fdata, function(adata) {
        
        # Nest data by county
        adata <- adata %>% nest(.by = county, .key = "ffdata")
        
        # Expand without bootstrap
        aadata <- adata %>% unnest(ffdata)
        
        # Run fixed-effect regression
        model <- feols(
          formula_object,
          aadata,
          weights = ~fraction,
          nthreads = 0,
          notes = FALSE,
          lean = TRUE
        )
        
        # Extract coefficients and confidence intervals
        tidy_res <- tidy(model)
        conf_int <- confint(model)
        tidy_res <- cbind(tidy_res, conf_int)
        
        # Extract goodness of fit (R²)
        glance_res <- glance(model)
        
        list(tidy_res = tidy_res, glance_res = glance_res)
      })
    )
}

# Run regression on full sample
full_sample_results <- run_full_sample_regression(
  data = fnl_data,
  formula_object = fml_base
)

# Extract coefficient table
coef_results <- full_sample_results %>%
  mutate(tidy_res = map(coefs, "tidy_res")) %>%
  select(crop_parent, tidy_res) %>%
  unnest(tidy_res)

# Extract goodness of fit table
goodness_of_fit <- full_sample_results %>%
  mutate(glance_res = map(coefs, "glance_res")) %>%
  select(crop_parent, glance_res) %>%
  unnest(glance_res)

# Save results
write.csv(
  coef_results,
  file.path(output_dir, "full_sample_coefficients_with_CI.csv"),
  row.names = FALSE
)

write.csv(
  goodness_of_fit,
  file.path(output_dir, "full_sample_goodness_of_fit.csv"),
  row.names = FALSE
)