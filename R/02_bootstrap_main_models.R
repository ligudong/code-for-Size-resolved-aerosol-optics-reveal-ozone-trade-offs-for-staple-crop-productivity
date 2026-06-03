source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

bootstrap_output_dir <- "data/derived/bootstrap"
dir.create(
  bootstrap_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

n_bootstrap <- 500

run_county_bootstrap <- function(data, formula_object, n = 500) {
  
  data %>%
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
            
            feols(
              formula_object,
              aadata,
              weights = ~fraction,
              nthreads = 0,
              notes = FALSE,
              lean = TRUE
            ) %>%
              tidy()
          },
          .id = "id",
          .progress = TRUE
        )
      })
    ) %>%
    select(-fdata)
}

data <- get_data()

# Bootstrap for the main model with interaction terms
bootstrap_base <- run_county_bootstrap(
  data = data,
  formula_object = fml_base,
  n = n_bootstrap
)

saveRDS(
  bootstrap_base,
  file.path(bootstrap_output_dir, "boots_base.rds")
)

# Bootstrap for the main model without interaction terms
bootstrap_base_no_interaction <- run_county_bootstrap(
  data = data,
  formula_object = fml_base_no_interaction,
  n = n_bootstrap
)

saveRDS(
  bootstrap_base_no_interaction,
  file.path(bootstrap_output_dir, "boots_base_no_interaction.rds")
)

# Prepare coefficient data for diagnostic plotting
coef_data_no_interaction <- bootstrap_base_no_interaction %>%
  unnest(coefs) %>%
  filter(term %in% c("Faod", "AOT40")) %>%
  select(crop_parent, term, estimate, id)

coef_boxplot_no_interaction <- coef_data_no_interaction %>%
  ggplot(aes(x = term, y = estimate, fill = term)) +
  geom_boxplot(
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2
  ) +
  facet_wrap(~ crop_parent, scales = "free") +
  scale_fill_manual(
    values = c(
      "Faod" = "#1f77b4",
      "AOT40" = "#ff7f0e"
    )
  ) +
  labs(
    title = "Bootstrap coefficient distributions",
    x = "Variable",
    y = "Estimate",
    fill = "Variable"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(coef_boxplot_no_interaction)

