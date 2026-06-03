##################fAOD#####################
source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

figure_dir <- "figures/supplementary/si22"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

data <- get_data() %>%
  drop_na() %>%
  filter(crop %in% c("Rice(LR)", "Rice(SR&ER)"))

model_names <- c("GOSIF", "CSIF", "RTSIF")

model_colors <- c(
  "GOSIF" = "#1f77b4",
  "CSIF"  = scales::alpha("#2ca02c", 0.35),
  "RTSIF" = scales::alpha("#d62728", 0.35)
)

model_fills <- c(
  "GOSIF" = scales::alpha("#1f77b4", 0.15),
  "CSIF"  = scales::alpha("#2ca02c", 0.06),
  "RTSIF" = scales::alpha("#d62728", 0.06)
)

plot_faod_sif_response <- function(adata) {
  
  data_modeling <- na.omit(adata)
  
  faod_range <- range(data_modeling$Faod, na.rm = TRUE)
  
  if (
    !is.finite(faod_range[1]) ||
    !is.finite(faod_range[2]) ||
    faod_range[1] == faod_range[2]
  ) {
    warning("Invalid fAOD range.")
    return(NULL)
  }
  
  models <- list(
    GOSIF = fml_irrigation_heterogeneity_gosif,
    CSIF  = fml_irrigation_heterogeneity_csif,
    RTSIF = fml_irrigation_heterogeneity_rtsif
  ) %>%
    map(
      ~ feols(
        .x,
        data = data_modeling,
        cluster = ~city,
        weights = ~fraction,
        nthreads = 0,
        lean = FALSE
      )
    )
  
  faod_hist <- seq(faod_range[1], faod_range[2], by = 0.01)
  faod_mid <- round(quantile(data_modeling$Faod, 0.5, na.rm = TRUE), 2)
  faod_offset <- (faod_range[2] - faod_range[1]) * 0.02
  
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
        mutate(Faod_x = rep(faod_range[2], 10))
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
  
  faod_pred %>%
    mutate(model = fct_relevel(model, model_names)) %>%
    ggplot(aes(x = Faod_x, y = fit, color = model, fill = model)) +
    geom_ribbon(
      aes(ymin = lwr, ymax = upr),
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 1.5) +
    geom_text_repel(
      aes(label = model),
      data = faod_label,
      show.legend = FALSE,
      min.segment.length = 0,
      box.padding = 0.8,
      direction = "y",
      nudge_x = 0.02,
      max.overlaps = 100,
      family = "Roboto Condensed",
      size = 5,
      seed = 2022
    ) +
    geom_xsidehistogram(
      aes(x = Faod),
      data = adata,
      bins = 80,
      fill = "#D3D3D3",
      color = "black",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = model_colors) +
    scale_fill_manual(values = model_fills) +
    ylab("Percentage Change in SIF") +
    xlab("fAOD") +
    labs(subtitle = "Response of GOSIF, CSIF and RTSIF to fAOD") +
    theme_half_open(18, font_family = "Roboto Condensed") +
    background_grid() +
    ggside(x.pos = "bottom", collapse = "x") +
    scale_xsidey_continuous(labels = NULL, breaks = NULL) +
    scale_x_continuous(
      limits = c(faod_range[1], faod_range[2] + faod_offset),
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme(
      ggside.panel.scale = 0.45,
      legend.title = element_blank(),
      legend.position = "right",
      panel.grid = element_blank()
    )
}

plot_si22 <- plot_faod_sif_response(data)

print(plot_si22)

ggsave(
  file.path(figure_dir, "figureSI22_faod_irrigation_sif_dataset_response_rice.tif"),
  plot_si22,
  width = 6,
  height = 4.5,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)



#####################AOT40#####################
source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

figure_dir <- "figures/supplementary/si22"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

data <- get_data() %>%
  drop_na() %>%
  filter(crop %in% c("Rice(LR)", "Rice(SR&ER)"))

model_names <- c("GOSIF", "CSIF", "RTSIF")

model_colors <- c(
  "GOSIF" = "#1f77b4",
  "CSIF"  = scales::alpha("#2ca02c", 0.35),
  "RTSIF" = scales::alpha("#d62728", 0.35)
)

model_fills <- c(
  "GOSIF" = scales::alpha("#1f77b4", 0.15),
  "CSIF"  = scales::alpha("#2ca02c", 0.06),
  "RTSIF" = scales::alpha("#d62728", 0.06)
)

plot_aot40_sif_response <- function(adata) {
  
  data_modeling <- na.omit(adata)
  aot40_range <- range(data_modeling$AOT40, na.rm = TRUE)
  
  if (
    !is.finite(aot40_range[1]) ||
    !is.finite(aot40_range[2]) ||
    aot40_range[1] == aot40_range[2]
  ) {
    warning("Invalid AOT40 range.")
    return(NULL)
  }
  
  models <- list(
    GOSIF = fml_irrigation_heterogeneity_gosif,
    CSIF  = fml_irrigation_heterogeneity_csif,
    RTSIF = fml_irrigation_heterogeneity_rtsif
  ) %>%
    map(
      ~ feols(
        .x,
        data = data_modeling,
        cluster = ~city,
        weights = ~fraction,
        nthreads = 0,
        lean = FALSE
      )
    )
  
  aot40_hist <- seq(aot40_range[1], aot40_range[2], by = 1)
  aot40_mid <- round(quantile(data_modeling$AOT40, 0.5, na.rm = TRUE), 2)
  aot40_offset <- (aot40_range[2] - aot40_range[1]) * 0.02
  
  aot40_pred <- map_dfr(
    models,
    function(amodel) {
      relpred(
        amodel,
        tibble(AOT40 = aot40_hist),
        tibble(AOT40 = rep(aot40_mid, length(aot40_hist)))
      ) %>%
        mutate(AOT40_x = aot40_hist)
    },
    .id = "model"
  ) %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 100))
  
  aot40_label <- map_dfr(
    models,
    function(amodel) {
      relpred(
        amodel,
        tibble(AOT40 = rep(aot40_range[2], 10)),
        tibble(AOT40 = rep(aot40_mid, 10))
      ) %>%
        mutate(AOT40_x = rep(aot40_range[2], 10))
    },
    .id = "model"
  ) %>%
    group_by(model) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(
      across(c(fit, lwr, upr), ~ expm1(.x) * 100),
      AOT40_x = AOT40_x + aot40_offset
    )
  
  aot40_pred %>%
    mutate(model = fct_relevel(model, model_names)) %>%
    ggplot(aes(x = AOT40_x, y = fit, color = model, fill = model)) +
    geom_ribbon(
      aes(ymin = lwr, ymax = upr),
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 1.5) +
    geom_text_repel(
      aes(label = model),
      data = aot40_label,
      show.legend = FALSE,
      min.segment.length = 0,
      box.padding = 0.8,
      direction = "y",
      nudge_x = 0.02,
      max.overlaps = 100,
      family = "Roboto Condensed",
      size = 5,
      seed = 2022
    ) +
    geom_xsidehistogram(
      aes(x = AOT40, y = after_stat(count * 3)),
      data = adata,
      bins = 80,
      fill = "#D3D3D3",
      color = "black",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = model_colors) +
    scale_fill_manual(values = model_fills) +
    ylab("Percentage Change in SIF") +
    xlab("AOT40") +
    labs(subtitle = "Response of GOSIF, CSIF and RTSIF to AOT40") +
    theme_half_open(18, font_family = "Roboto Condensed") +
    background_grid() +
    ggside(x.pos = "bottom", collapse = "x") +
    scale_xsidey_continuous(labels = NULL, breaks = NULL) +
    scale_x_continuous(
      limits = c(aot40_range[1], aot40_range[2] + aot40_offset),
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme(
      ggside.panel.scale = 3,
      legend.title = element_blank(),
      legend.position = "right",
      panel.grid = element_blank()
    )
}

plot_si22_aot40 <- plot_aot40_sif_response(data)

print(plot_si22_aot40)

ggsave(
  file.path(figure_dir, "figureSI22_aot40_irrigation_sif_dataset_response_rice.tif"),
  plot_si22_aot40,
  width = 6,
  height = 4.5,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)