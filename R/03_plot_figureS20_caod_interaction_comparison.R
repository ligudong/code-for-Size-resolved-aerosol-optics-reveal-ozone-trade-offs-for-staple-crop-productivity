# Load necessary packages and functions
source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")  # Contains fml_base4 and fml_base_no_interaction4

# Create output directory
figure_dir <- "figures/supplementary/si20"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Load and preprocess data
data <- get_data() %>%
  drop_na()

# Filter for combined Rice crops
data <- data %>% filter(crop %in% c("Rice(LR)", "Rice(SR&ER)"))

# Define model names
model_names <- c("interaction", "no_interaction")

# Define the plotting function
get_plot <- function(adata) {
  message("Processing crop data")
  message("Caod summary: ", summary(adata$Caod))
  
  data_modeling <- na.omit(adata)
  Caod_range <- range(data_modeling$Caod, na.rm = TRUE)
  
  if (!is.finite(Caod_range[1]) || !is.finite(Caod_range[2]) || Caod_range[1] == Caod_range[2]) {
    warning("Invalid Caod range, skip plot")
    return(NULL)
  }
  
  # Fit two models: interaction and no interaction
  models <- lst(
    fml_base4, 
    fml_base_no_interaction4
  ) %>%
    map(
      ~ feols(
        .x,
        data = data_modeling,
        cluster = ~ city,
        weights = ~ fraction,
        nthreads = 0,
        lean = TRUE
      ),
      .progress = TRUE
    ) %>%
    set_names(model_names)
  
  # Construct Caod response axis
  Caod_hist <- seq(Caod_range[1], Caod_range[2], 0.001)
  Caod_mid <- round(quantile(data_modeling$Caod, 0.5, na.rm = TRUE), 2)
  
  # Generate predicted values
  Caod_pred <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(Caod = Caod_hist),
      tibble(Caod = rep(Caod_mid, length(Caod_hist)))
    ) %>%
      mutate(Caod_x = Caod_hist)
  }, .id = "model") %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 100))
  
  # Offset for label positioning
  Caod_offset <- (Caod_range[2] - Caod_range[1]) * 0.01
  
  # Label positions
  Caod_lab <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(Caod = rep(Caod_range[2], 10)),
      tibble(Caod = rep(Caod_mid, 10))
    ) %>%
      mutate(Caod_x = rep(Caod_range[2], 10))
  }, .id = "model") %>%
    group_by(model) %>% slice_head(n = 1) %>% ungroup() %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 100),
           Caod_x = Caod_x + Caod_offset)
  
  # Define colors and fills
  model_colors <- c(
    "interaction" = "#1f77b4",
    "no_interaction" = scales::alpha("#ff7f0e", 0.40)
  )
  model_fills <- c(
    "interaction" = scales::alpha("#1f77b4", 0.15),
    "no_interaction" = scales::alpha("#ff7f0e", 0.06)
  )
  
  # Plotting
  p <- Caod_pred %>%
    mutate(model = fct_relevel(model, model_names)) %>%
    ggplot(aes(x = Caod_x, y = fit, color = model, fill = model)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), color = NA, show.legend = FALSE) +
    geom_line(size = 1.5) +
    geom_text_repel(aes(label = model),
                    data = Caod_lab,
                    show.legend = FALSE,
                    min.segment.length = 0,
                    box.padding = 1,
                    direction = "y",
                    nudge_x = 0.02,
                    nudge_y = 0.02,
                    max.overlaps = 100,
                    family = "Roboto Condensed",
                    size = 5,
                    seed = 2022) +
    geom_xsidehistogram(aes(x = Caod), data = adata, bins = 80,
                        fill = "#D3D3D3", color = "black", size = 0.2, inherit.aes = FALSE) +
    scale_color_manual(values = model_colors) +
    scale_fill_manual(values = model_fills) +
    ylab("Percentage Change in SIF") +
    xlab("cAOD") +
    labs(subtitle = "Crop Response Plot") +
    theme_half_open(18, font_family = "Roboto Condensed") +
    background_grid() +
    ggside(x.pos = "bottom", collapse = "x") +
    scale_xsidey_continuous(labels = NULL, breaks = NULL) +
    scale_x_continuous(limits = c(Caod_range[1], Caod_range[2] + Caod_offset),
                       expand = expansion(mult = c(0, 0.15))) +
    theme(ggside.panel.scale = 0.5,
          legend.title = element_blank(),
          legend.position = "right",
          panel.grid = element_blank())
  
  return(p)
}

# Generate the plot
plots <- get_plot(data)

# Create figure directory if not exists
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

# Print and save the plot
print(plots)
ggsave(
  filename = file.path(figure_dir, "figureSI20_caod_response_Rice.tif"),
  plot = plots,
  width = 6,
  height = 4.5,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)