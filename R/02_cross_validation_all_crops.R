source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

data_dir   <- "data/derived/crop/maize"
output_dir <- "data/derived/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fnl_data <- get_data()
n_folds <- 10  # Number of folds

crops <- c("Maize", "Rice", "Wheat")
formulas <- list(
  base = fml_base,
  base_no_interaction = fml_base_no_interaction
)

# Function to perform cross-validation per crop and per formula
run_cross_validation <- function(data, formula_object, n_folds = 10) {
  
  data %>%
    nest(fdata = -crop_parent) %>%
    mutate(coefs = map(fdata, function(adata) {
      
      adata <- adata %>% nest(.by = county, .key = "ffdata")
      
      map_dfr(1:n_folds, function(fold) {
        
        set.seed(fold)
        
        aadata <- adata %>%
          slice_sample(prop = 1, replace = TRUE) %>%
          unnest(ffdata)
        
        # Demean to remove city/year fixed effects
        aadata <- aadata %>%
          group_by(city, year) %>%
          mutate(
            across(
              c(Caod, Faod, cloud, AOT40, bin1:bin42, root_1:root_9, d1:d13),
              ~ . - mean(.),
              .names = "demeaned_{.col}"
            )
          ) %>%
          ungroup()
        
        model <- feols(
          formula_object,
          aadata,
          weights = ~fraction,
          nthreads = 0,
          notes = FALSE,
          lean = TRUE,
          cluster = ~x_y
        )
        
        tidy_res <- tidy(model)
        conf_int <- confint(model)
        tidy_res <- cbind(tidy_res, conf_int)
        glance_res <- glance(model)
        
        tidy_res %>%
          mutate(
            r_squared = glance_res$r.squared,
            adj_r_squared = glance_res$adj.r.squared,
            within_r_squared = glance_res$within.r.squared,
            fold_id = fold
          )
        
      }, .id = "id", .progress = TRUE)
    }))
}

# Run cross-validation for each formula and crop
for (formula_name in names(formulas)) {
  
  formula_obj <- formulas[[formula_name]]
  
  cv_results <- map_dfr(crops, function(crop_name) {
    crop_data <- fnl_data %>% filter(crop == crop_name)
    run_cross_validation(crop_data, formula_obj, n_folds)
  })
  
  # Save results
  saveRDS(
    cv_results,
    file.path(output_dir, paste0("cross_validation_", formula_name, ".rds"))
  )
}