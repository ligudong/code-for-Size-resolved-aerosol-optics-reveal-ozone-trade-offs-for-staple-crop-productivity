# Load packages and helper functions
source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")  # contains fml_base and fml_base_no_interaction

library(dplyr)
library(ggplot2)
library(fixest)
library(forcats)
library(purrr)
library(ggrepel)
library(ggside)
library(tibble)

# Output directory
figure_dir <- "figures/supplementary/si8"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
data <- get_data() %>% drop_na()

# Filter per crop
data_maize <- data %>% filter(crop == "Maize")
data_rice  <- data %>% filter(crop %in% c("Rice(LR)","Rice(SR&ER)"))
data_wheat <- data %>% filter(crop == "Wheat")

# Define model names
model_names <- c("base", "no_interaction")

# Function to generate AOT40 response plot
get_aot40_response_plot <- function(adata) {
  
  message("Processing crop data")
  message("AOT40 summary: ", summary(adata$AOT40))
  
  data_modeling <- na.omit(adata)
  
  # Check range
  AOT40_range <- range(data_modeling$AOT40, na.rm = TRUE)
  if (!is.finite(AOT40_range[1]) || !is.finite(AOT40_range[2]) || AOT40_range[1]==AOT40_range[2]) {
    warning("Invalid AOT40 range, skip plot")
    return(NULL)
  }
  
  # Fit models
  models <- lst(
    fml_base,
    fml_base_no_interaction
  ) %>%
    map(~ feols(.x, data_modeling, cluster = ~city, weights = ~fraction, nthreads = 0, lean = TRUE),
        .progress = TRUE
    ) %>%
    set_names(model_names)
  
  # Prediction axis
  AOT40_hist <- seq(AOT40_range[1], AOT40_range[2], 1)
  AOT40_mid  <- round(quantile(adata$AOT40, 0.5, na.rm=TRUE), 2)
  AOT40_offset <- (AOT40_range[2] - AOT40_range[1]) * 0.01
  
  # Calculate predicted values
  AOT40_df <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(AOT40 = AOT40_hist),
      tibble(AOT40 = rep(AOT40_mid, length(AOT40_hist)))
    ) %>%
      mutate(AOT40_x = AOT40_hist)
  }, .id = "model") %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x)*1e2))
  
  # Labels for plotting
  AOT40_lab <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(AOT40 = rep(AOT40_range[2], 10)),
      tibble(AOT40 = rep(AOT40_mid, 10))
    ) %>%
      mutate(AOT40_x = rep(AOT40_range[2], 10))
  }, .id = "model") %>%
    group_by(model) %>%
    slice_head(n=1) %>%
    ungroup() %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x)*1e2)) %>%
    mutate(AOT40_x = AOT40_x + AOT40_offset)
  
  # Plot
  model_colors <- c("base"="#1f77b4","no_interaction"="#ff7f0e")
  
  p1 <- AOT40_df %>%
    mutate(model = fct_relevel(model, rev(model_names))) %>%
    ggplot(aes(x=AOT40_x, y=fit, color=model)) +
    geom_ribbon(aes(x=AOT40_x, ymin=lwr, ymax=upr, fill=model),
                inherit.aes=FALSE, alpha=0.2, show.legend=FALSE,
                size=0.01, data=. %>% filter(model=="base")) +
    geom_line(data=. %>% filter(model=="base"), size=1.5) +
    geom_line(data=. %>% filter(model=="no_interaction"), size=1.5) +
    scale_color_manual(values=model_colors) +
    geom_text_repel(
      data=AOT40_lab,
      aes(label=model),
      min.segment.length=0,
      box.padding=1,
      direction="y",
      nudge_x=0.02,
      nudge_y=0.02,
      xlim=c(AOT40_range[2]+AOT40_offset, AOT40_range[2]*1.1),
      max.overlaps=100,
      family="Roboto Condensed",
      size=5,
      seed=2022
    ) +
    ylab("Percentage Change in SIF") +
    xlab("AOT40") +
    theme_half_open(18, font_family="Roboto Condensed") +
    background_grid() +
    theme(
      ggside.panel.scale=0.2,
      legend.title=element_blank(),
      legend.position=c(1.2,0.5),
      panel.grid=element_blank()
    ) +
    geom_xsidehistogram(aes(x=AOT40),
                        data=adata,
                        bins=80,
                        fill="#D3D3D3",
                        color="black",
                        size=0.2,
                        inherit.aes=FALSE) +
    ggside(x.pos="bottom", collapse="x") +
    scale_xsidey_continuous(labels=NULL, breaks=NULL) +
    scale_x_continuous(limits=c(AOT40_range[1], AOT40_range[2]+AOT40_offset),
                       expand=expansion(mult=c(0,0.15))) +
    labs(subtitle="Crop Response to AOT40")
  
  return(p1)
}

# Generate plots per crop
plot_maize <- get_aot40_response_plot(data %>% filter(crop=="Maize"))
plot_rice  <- get_aot40_response_plot(data %>% filter(crop %in% c("Rice(LR)","Rice(SR&ER)")))
plot_wheat <- get_aot40_response_plot(data %>% filter(crop=="Wheat"))

# Save TIFF
ggsave(file.path(figure_dir,"si8_aot40_response_maize.tif"), plot_maize, width=6, height=4.5, dpi=300, device="tiff")
ggsave(file.path(figure_dir,"si8_aot40_response_rice.tif"),  plot_rice,  width=6, height=4.5, dpi=300, device="tiff")
ggsave(file.path(figure_dir,"si8_aot40_response_wheat.tif"), plot_wheat, width=6, height=4.5, dpi=300, device="tiff")