source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

figure_dir <- "figures/supplementary/si8"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

data <- get_data() %>%
  drop_na() %>%
  filter(crop == "Maize")

model_names <- c("base", "no_interaction")

get_faod_response_plot <- function(adata) {
  
  data_modeling <- na.omit(adata)
  
  faod_range <- range(data_modeling$Faod, na.rm = TRUE)
  
  if (!is.finite(faod_range[1]) ||
      !is.finite(faod_range[2]) ||
      faod_range[1] == faod_range[2]) {
    warning("Invalid fAOD range.")
    return(NULL)
  }
  
  models <- list(
    base = fml_base,
    no_interaction = fml_base_no_interaction
  ) %>%
    map(
      ~ feols(
        .x,
        data_modeling,
        cluster = ~ city,
        weights = ~ fraction,
        nthreads = 0,
        lean = TRUE
      )
    )
  
  faod_hist <- seq(faod_range[1], faod_range[2], by = 0.01)
  faod_mid <- round(quantile(data_modeling$Faod, 0.5, na.rm = TRUE), 2)
  faod_offset <- (faod_range[2] - faod_range[1]) * 0.01
  
  faod_pred <- map_dfr(
    models,
    function(amodel) {
      relpred(
        amodel,
        tibble(Faod = faod_hist),
        tibble(Faod = rep(faod_mid, length(faod_hist)))
      ) %>%
        mutate(Faod_x = faod_hist)
    },
    .id = "model"
  ) %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 100))
  
  faod_label <- map_dfr(
    models,
    function(amodel) {
      relpred(
        amodel,
        tibble(Faod = rep(faod_range[2], 10)),
        tibble(Faod = rep(faod_mid, 10))
      ) %>%
        mutate(Faod_x = faod_range[2])
    },
    .id = "model"
  ) %>%
    group_by(model) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(
      across(c(fit, lwr, upr), ~ expm1(.x) * 100),
      Faod_x = Faod_x + faod_offset
    )
  
  model_colors <- c(
    "base" = "#1f77b4",
    "no_interaction" = "#ff7f0e"
  )
  
  faod_pred %>%
    mutate(model = fct_relevel(model, rev(model_names))) %>%
    ggplot(aes(x = Faod_x, y = fit, color = model)) +
    geom_ribbon(
      data = . %>% filter(model == "base"),
      aes(x = Faod_x, ymin = lwr, ymax = upr, fill = model),
      inherit.aes = FALSE,
      alpha = 0.2,
      linewidth = 0.01,
      show.legend = FALSE
    ) +
    geom_line(
      data = . %>% filter(model == "base"),
      linewidth = 1.5,
      linetype = "solid"
    ) +
    geom_line(
      data = . %>% filter(model == "no_interaction"),
      linewidth = 1.5,
      linetype = "solid"
    ) +
    scale_color_manual(values = model_colors) +
    scale_fill_manual(values = model_colors) +
    geom_text_repel(
      data = faod_label,
      aes(label = model),
      show.legend = FALSE,
      min.segment.length = 0,
      box.padding = 1,
      direction = "y",
      nudge_x = 0.02,
      nudge_y = 0.02,
      max.overlaps = 100,
      family = "Roboto Condensed",
      size = 5,
      seed = 2022
    ) +
    geom_xsidehistogram(
      data = adata,
      aes(x = Faod),
      bins = 80,
      fill = "#D3D3D3",
      color = "black",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    ggside(x.pos = "bottom", collapse = "x") +
    scale_x_continuous(
      limits = c(faod_range[1], faod_range[2] + faod_offset),
      expand = expansion(mult = c(0, 0.15))
    ) +
    scale_xsidey_continuous(labels = NULL, breaks = NULL) +
    ylab("Percentage change in SIF") +
    xlab("fAOD") +
    theme_half_open(18, font_family = "Roboto Condensed") +
    background_grid() +
    theme(
      ggside.panel.scale = 0.5,
      legend.title = element_blank(),
      legend.position = c(1.2, 0.5),
      panel.grid = element_blank()
    )
}

plot_faod <- get_faod_response_plot(data)

print(plot_faod)

ggsave(
  file.path(figure_dir, "si8_faod_response_maize.tif"),
  plot_faod,
  width = 6,
  height = 4.5,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(figure_dir, "si8_faod_response_maize.pdf"),
  plot_faod,
  width = 6,
  height = 4.5
)