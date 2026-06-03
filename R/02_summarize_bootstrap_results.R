source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

bootstrap_output_dir <- "data/derived/bootstrap"
summary_output_dir <- "data/derived/results"

dir.create(summary_output_dir, recursive = TRUE, showWarnings = FALSE)

n_bootstrap <- 500

run_bootstrap_with_summary <- function(data, formula_object, n = 500) {
  
  bootstrap_results <- data %>%
    nest(fdata = -crop_parent) %>%
    mutate(
      coefs = map(fdata, function(adata) {
        
        adata <- adata %>%
          nest(.by = county, .key = "ffdata")
        
        map_dfr(
          1:n,
          function(anum) {
            
            set.seed(anum)
            
            aadata <- adata %>%
              slice_sample(prop = 1, replace = TRUE) %>%
              unnest(ffdata)
            
            model <- feols(
              formula_object,
              aadata,
              weights = ~fraction,
              nthreads = 0,
              notes = FALSE,
              lean = TRUE,
              cluster = ~x_y
            )
            
            tidy(model) %>%
              mutate(
                r_squared = glance(model)$r.squared
              )
          },
          .id = "id",
          .progress = TRUE
        )
      })
    ) %>%
    select(-fdata)
  
  bootstrap_summary <- bootstrap_results %>%
    mutate(
      coefs_summary = map(coefs, function(coef_data) {
        coef_data %>%
          group_by(term) %>%
          summarise(
            mean_estimate = mean(estimate, na.rm = TRUE),
            sd_estimate = sd(estimate, na.rm = TRUE),
            ci_low = quantile(estimate, 0.025, na.rm = TRUE),
            ci_high = quantile(estimate, 0.975, na.rm = TRUE),
            mean_r_squared = mean(r_squared, na.rm = TRUE),
            .groups = "drop"
          )
      })
    ) %>%
    select(crop_parent, coefs_summary) %>%
    unnest(coefs_summary)
  
  list(
    bootstrap_results = bootstrap_results,
    bootstrap_summary = bootstrap_summary
  )
}

data <- get_data()

base_output <- run_bootstrap_with_summary(
  data = data,
  formula_object = fml_base,
  n = n_bootstrap
)

base_no_interaction_output <- run_bootstrap_with_summary(
  data = data,
  formula_object = fml_base_no_interaction,
  n = n_bootstrap
)

saveRDS(
  base_output$bootstrap_results,
  file.path(bootstrap_output_dir, "boots_base_clustered.rds")
)

saveRDS(
  base_no_interaction_output$bootstrap_results,
  file.path(bootstrap_output_dir, "boots_base_no_interaction_clustered.rds")
)

saveRDS(
  base_output$bootstrap_summary,
  file.path(summary_output_dir, "bootstrap_summary_base.rds")
)

saveRDS(
  base_no_interaction_output$bootstrap_summary,
  file.path(summary_output_dir, "bootstrap_summary_base_no_interaction.rds")
)