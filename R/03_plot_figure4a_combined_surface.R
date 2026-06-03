source("R/00_load_packages.R")
source("R/00_helper_functions.R")

impact_dir <- "data/derived/impacts"
output_dir <- "figures/figure4"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

faod_impact_file <- file.path(impact_dir, "impacts_Faod.qs")
ozone_impact_file <- file.path(impact_dir, "impacts_ozone_AOT40.qs")

faod_impacts <- qread(faod_impact_file, nthreads = qn)
ozone_impacts <- qread(ozone_impact_file, nthreads = qn)

extract_faod_response <- function(data, crop_name) {
  data %>%
    filter(crop_parent == crop_name) %>%
    select(crop_parent, year, Faod_results) %>%
    unnest(Faod_results) %>%
    group_by(crop_parent, faod_level) %>%
    summarise(
      value = mean(`50%`, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    rename(Faod_level = faod_level)
}

extract_ozone_response <- function(data, crop_name) {
  data %>%
    filter(crop_parent == crop_name) %>%
    select(crop_parent, year, ozone_results) %>%
    unnest(ozone_results) %>%
    group_by(crop_parent, peak_level) %>%
    summarise(
      value = mean(`50%`, na.rm = TRUE) * 100,
      .groups = "drop"
    )
}

plot_combined_surface <- function(crop_name, output_file, zlim = c(-50, 30)) {
  
  faod_response <- extract_faod_response(faod_impacts, crop_name)
  ozone_response <- extract_ozone_response(ozone_impacts, crop_name)
  
  combined_data <- expand.grid(
    Faod_level = faod_response$Faod_level,
    peak_level = ozone_response$peak_level
  ) %>%
    left_join(faod_response, by = "Faod_level") %>%
    rename(value_faod = value) %>%
    left_join(
      ozone_response %>% select(peak_level, value),
      by = "peak_level"
    ) %>%
    rename(value_ozone = value) %>%
    mutate(
      combined_value = value_faod + value_ozone
    )
  
  colors_blue_to_red <- c(
    "#061178", "#10239E", "#1D39C4", "#4465EB", "#6682F5",
    "#8BA2FF", "#ADC6FF", "#D6E4FF", "#F0F5FF",
    "#FFE4D6", "#FFD2B9", "#FFB186", "#FF8F50", "#FF762A",
    "#F55C08", "#D94F03", "#AD2102", "#871400"
  )
  
  color_positions <- rescale(
    c(
      seq(zlim[1], 0, length.out = 9),
      seq(1, zlim[2], length.out = 9)
    ),
    from = zlim
  )
  
  p <- ggplot(
    combined_data,
    aes(
      x = Faod_level,
      y = peak_level,
      z = combined_value
    )
  ) +
    geom_tile(aes(fill = combined_value)) +
    geom_contour(color = "grey70", linewidth = 0.4) +
    geom_contour(breaks = 0, color = "black", linewidth = 0.8) +
    scale_fill_gradientn(
      colors = colors_blue_to_red,
      values = color_positions,
      limits = zlim,
      oob = scales::squish,
      name = "SIF change (%)"
    ) +
    labs(
      x = "Annual fAOD level",
      y = expression(Annual~O[3]~level),
      title = crop_name
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
  
  ggsave(
    output_file,
    p,
    width = 6,
    height = 5,
    dpi = 600,
    compression = "lzw"
  )
  
  p
}

plot_combined_surface(
  crop_name = "Maize",
  output_file = file.path(output_dir, "figure4a_faod_ozone_surface_maize.tif"),
  zlim = c(-50, 30)
)

plot_combined_surface(
  crop_name = "Rice",
  output_file = file.path(output_dir, "figure4a_faod_ozone_surface_rice.tif"),
  zlim = c(-50, 30)
)

plot_combined_surface(
  crop_name = "Wheat",
  output_file = file.path(output_dir, "figure4a_faod_ozone_surface_wheat.tif"),
  zlim = c(-50, 30)
)