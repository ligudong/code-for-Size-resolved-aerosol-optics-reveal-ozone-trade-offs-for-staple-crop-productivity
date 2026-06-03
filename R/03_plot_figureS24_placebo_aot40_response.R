source("scripts/00_load_packages.R")
source("scripts/loadFunctions.R")
source("scripts/loadFormulas.R") 

# Output directory
figure_dir <- "figures/supplementary/si24_si25"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Load and preprocess data
data_all <- get_data() %>%
  drop_na()

# Placebo test parameters
set.seed(123)
n_simulations <- 500

# Crop list
crops <- c("Rice(LR)", "Rice(SR&ER)", "Maize", "Wheat")

# Loop through crops for placebo tests
for(crop_name in crops){
  
  adata <- data_all %>% filter(crop %in% crop_name)
  
  # Initialize result container
  placebo_results <- data.frame(beta = numeric(0))
  count_valid <- 0
  
  while(count_valid < n_simulations){
    # Shuffle AOT40 and Faod for interaction test
    adata_placebo <- adata %>%
      mutate(
        Faod_placebo = sample(Faod),
        AOT40_placebo = sample(AOT40),
        Faod_AOT40_placebo = Faod_placebo * AOT40_placebo
      ) %>%
      drop_na(GOSIF_sum)  # Ensure dependent variable exists
    
    # Skip if insufficient rows
    if(nrow(adata_placebo) < 10) next
    
    # Construct regression formula
    fml_placebo <- str_c(
      "GOSIF_sum ~ ",
      c(
        str_c(c(
          "Faod_placebo", "Faod_placebo^2",
          "Caod", "Caod^2",
          "AOT40_placebo", "cloud", "cloud^2",
          "co2",
          str_c("bin", 1:42),
          str_c("root_", 1:9),
          str_c("d", 1:13)
        ), "* irg_fraction * VPD")
      ) %>% str_flatten(collapse = "+"),
      "| x_y[year]"
    ) %>% as.formula()
    
    # Run regression
    model_placebo <- tryCatch(
      feols(fml_placebo, data = adata_placebo, weights = ~fraction, notes = FALSE, lean = TRUE),
      error = function(e) NULL
    )
    
    # Store coefficient if model valid
    if(!is.null(model_placebo)){
      coef_val <- coef(model_placebo)["Faod_AOT40_placebo"]
      if(!is.na(coef_val)){
        placebo_results <- rbind(placebo_results, data.frame(beta = coef_val))
        count_valid <- count_valid + 1
      }
    }
  }
  
  # Create placebo histogram plot
  p_placebo <- ggplot(placebo_results, aes(x = beta)) +
    geom_histogram(binwidth = 0.000005, fill = "lightblue", color = "darkblue", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    labs(x = "Placebo FaOD*AOT40 Coefficient", y = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "plain"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(limits = c(-0.0001, 0.0001), breaks = seq(-0.0001, 0.0001, 0.00005), expand = c(0,0))
  
  # Print and save plot
  print(p_placebo)
  
  ggsave(
    filename = file.path(figure_dir, paste0("placebo_FaodAOT40_bothshuffled_", gsub("\\(|\\)", "", crop_name), ".tif")),
    plot = p_placebo,
    width = 6,
    height = 4.5,
    dpi = 300,
    device = "tiff"
  )
}