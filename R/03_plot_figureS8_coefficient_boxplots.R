source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

figure_dir <- "figures/supplementary/si8"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

data <- get_data() %>%
  drop_na()

n_bootstrap <- 500

models <- list(
  Base = fml_base,
  No_interaction = fml_base_no_interaction,
  Base_linear = fml_base_lin,
  No_ozone = fml_no_O3
)

target_terms <- c("Faod", "Caod", "AOT40")

calculate_coefficients <- function(data, formula_object, n_iter, model_name, target_term) {
  
  data %>%
    nest(fdata = -crop_parent) %>%
    mutate(
      coefs = map(fdata, function(adata) {
        
        adata <- adata %>%
          nest(.by = county, .key = "ffdata")
        
        map_dfr(
          1:n_iter,
          function(anum) {
            
            set.seed(anum)
            
            aadata <- adata %>%
              slice_sample(prop = 1, replace = TRUE) %>%
              unnest(ffdata)
            
            feols(
              formula_object,
              aadata,
              weights = ~ fraction,
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
    unnest(coefs) %>%
    filter(term == target_term) %>%
    select(crop_parent, term, estimate, id) %>%
    mutate(model = model_name)
}

make_boxplot_data <- function(target_term) {
  
  imap_dfr(
    models,
    function(formula_object, model_name) {
      calculate_coefficients(
        data = data,
        formula_object = formula_object,
        n_iter = n_bootstrap,
        model_name = model_name,
        target_term = target_term
      )
    }
  )
}

plot_coefficient_boxplot <- function(coef_data, target_term) {
  
  coef_data %>%
    mutate(
      model = factor(
        model,
        levels = c("Base", "No_interaction", "Base_linear", "No_ozone")
      )
    ) %>%
    ggplot(aes(x = model, y = estimate, fill = model)) +
    geom_boxplot(
      outlier.color = "red",
      outlier.shape = 16,
      outlier.size = 1.5,
      alpha = 0.65
    ) +
    geom_jitter(
      aes(color = model),
      width = 0.18,
      size = 1.1,
      alpha = 0.35
    ) +
    facet_wrap(~ crop_parent, scales = "free_y") +
    scale_fill_manual(
      values = c(
        "Base" = "#1f77b4",
        "No_interaction" = "#ff7f0e",
        "Base_linear" = "#2ca02c",
        "No_ozone" = "#d62728"
      )
    ) +
    scale_color_manual(
      values = c(
        "Base" = "#1f77b4",
        "No_interaction" = "#ff7f0e",
        "Base_linear" = "#2ca02c",
        "No_ozone" = "#d62728"
      )
    ) +
    labs(
      title = paste0(target_term, " coefficient distribution across models"),
      x = NULL,
      y = "Estimate",
      fill = NULL,
      color = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
}

coef_results <- map(
  target_terms,
  function(target_term) {
    
    coef_data <- make_boxplot_data(target_term)
    
    saveRDS(
      coef_data,
      file.path(
        figure_dir,
        paste0("si8_", tolower(target_term), "_coefficient_bootstrap.rds")
      )
    )
    
    p <- plot_coefficient_boxplot(coef_data, target_term)
    
    ggsave(
      file.path(
        figure_dir,
        paste0("si8_", tolower(target_term), "_coefficient_boxplot.tif")
      ),
      p,
      width = 12,
      height = 8,
      dpi = 600,
      device = "tiff",
      compression = "lzw"
    )
    
    p
  }
)

names(coef_results) <- target_terms