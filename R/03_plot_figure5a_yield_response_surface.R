source("R/00_load_packages.R")
source("R/00_helper_functions.R")

library(dplyr)
library(ggplot2)
library(scales)
library(latex2exp)
library(ggstar)
library(qs)

input_file <- "data/derived/yield_response/yield_response_surface.rds"
output_dir <- "figures/yield_response"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

yield_response <- readRDS(input_file)

plot_yield_heatmap <- function(data, crop_name, output_name) {
  
  df_filtered <- data %>%
    filter(
      faod_level >= 0,
      faod_level <= 1.5,
      peak_level >= 30,
      peak_level <= 170
    )
  
  p <- ggplot(
    df_filtered,
    aes(
      x = faod_level,
      y = peak_level,
      fill = value
    )
  ) +
    geom_tile() +
    scale_fill_gradientn(
      name = TeX("Percentage change in yield"),
      limits = c(-5, 5),
      oob = squish,
      colours = c(
        "#d73027",
        "#fdae61",
        "white",
        "#d9ef8b",
        "#1a9850"
      )
    ) +
    geom_contour(
      aes(z = value),
      color = "grey",
      bins = 20
    ) +
    scale_x_continuous(limits = c(0, 1.5)) +
    scale_y_continuous(limits = c(30, 170)) +
    labs(
      title = crop_name,
      x = "fAOD level",
      y = "Ozone peak level"
    ) +
    geom_point(
      data = df_filtered %>%
        filter(faod_level == 0.4, peak_level == 100),
      aes(x = faod_level, y = peak_level),
      color = "black",
      size = 3,
      shape = 16
    ) +
    geom_star(
      data = df_filtered %>%
        slice_max(value, n = 1, with_ties = FALSE),
      aes(x = faod_level, y = peak_level),
      fill = "yellow",
      color = NA,
      size = 3
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(output_dir, output_name),
    plot = p,
    width = 8,
    height = 7,
    dpi = 300,
    units = "in",
    device = "tiff"
  )
  
  p
}

plot_yield_heatmap(
  data = yield_response %>% filter(crop == "Maize"),
  crop_name = "Maize",
  output_name = "heatmap_yield_response_maize.tif"
)

plot_yield_heatmap(
  data = yield_response %>% filter(crop == "Rice"),
  crop_name = "Rice",
  output_name = "heatmap_yield_response_rice.tif"
)

plot_yield_heatmap(
  data = yield_response %>% filter(crop == "Wheat"),
  crop_name = "Wheat",
  output_name = "heatmap_yield_response_wheat.tif"
)

plot_yield_heatmap(
  data = yield_response %>%
    group_by(faod_level, peak_level) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop"),
  crop_name = "All crops",
  output_name = "heatmap_yield_response_all_crops.tif"
)