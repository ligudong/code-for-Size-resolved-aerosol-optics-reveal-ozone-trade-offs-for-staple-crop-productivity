source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

figure_dir <- "figures/supplementary/si23"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

data_all <- get_data() %>%
  drop_na()

model_colors <- c(
  "Faod × AOT40" = "#1f77b4",
  "No Faod × AOT40" = scales::alpha("#ff7f0e", 0.40)
)

model_fills <- c(
  "Faod × AOT40" = scales::alpha("#1f77b4", 0.15),
  "No Faod × AOT40" = scales::alpha("#ff7f0e", 0.06)
)

make_newdata <- function(var, value) {
  as_tibble(setNames(list(value), var))
}

plot_interaction_response <- function(
    adata,
    var,
    x_label,
    subtitle_text,
    hist_step,
    hist_scale = 1
) {
  
  data_modeling <- na.omit(adata)
  var_range <- range(data_modeling[[var]], na.rm = TRUE)
  
  if (
    !is.finite(var_range[1]) ||
    !is.finite(var_range[2]) ||
    var_range[1] == var_range[2]
  ) {
    warning(paste("Invalid range for", var))
    return(NULL)
  }
  
  models <- list(
    `Faod × AOT40` = fml_irrigation_heterogeneity_gosif,
    `No Faod × AOT40` = fml_irrigation_heterogeneity_gosif_no_interaction
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
  
  x_hist <- seq(var_range[1], var_range[2], by = hist_step)
  x_mid <- round(quantile(data_modeling[[var]], 0.5, na.rm = TRUE), 2)
  x_offset <- (var_range[2] - var_range[1]) * 0.02
  
  pred_df <- map_dfr(
    models,
    function(amodel) {
      relpred(
        amodel,
        make_newdata(var, x_hist),
        make_newdata(var, rep(x_mid, length(x_hist)))
      ) %>%
        mutate(x_value = x_hist)
    },
    .id = "model"
  ) %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x) * 100))
  
  label_df <- map_dfr(
    models,
    function(amodel) {
      relpred(
        amodel,
        make_newdata(var, rep(var_range[2], 10)),
        make_newdata(var, rep(x_mid, 10))
      ) %>%
        mutate(x_value = rep(var_range[2], 10))
    },
    .id = "model"
  ) %>%
    group_by(model) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(
      across(c(fit, lwr, upr), ~ expm1(.x) * 100),
      x_value = x_value + x_offset
    )
  
  ggplot(pred_df, aes(x = x_value, y = fit, color = model, fill = model)) +
    geom_ribbon(
      aes(ymin = lwr, ymax = upr),
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 1.5) +
    geom_text_repel(
      aes(label = model),
      data = label_df,
      show.legend = FALSE,
      box.padding = 0.8,
      direction = "y",
      nudge_x = 0.02,
      family = "Roboto Condensed",
      size = 5,
      seed = 2022
    ) +
    geom_xsidehistogram(
      aes(x = .data[[var]], y = after_stat(count * hist_scale)),
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
    xlab(x_label) +
    labs(subtitle = subtitle_text) +
    theme_half_open(18, font_family = "Roboto Condensed") +
    background_grid() +
    ggside(x.pos = "bottom", collapse = "x") +
    scale_xsidey_continuous(labels = NULL, breaks = NULL) +
    scale_x_continuous(
      limits = c(var_range[1], var_range[2] + x_offset),
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme(
      ggside.panel.scale = ifelse(var == "AOT40", 3, 0.45),
      legend.title = element_blank(),
      legend.position = "right",
      panel.grid = element_blank()
    )
}

plot_specs <- tribble(
  ~crop_name, ~crop_filter, ~var, ~x_label, ~hist_step, ~hist_scale,
  "rice",  list(c("Rice(LR)", "Rice(SR&ER)")), "Faod",  "fAOD",  0.01, 1,
  "maize", list("Maize"),                        "AOT40", "AOT40", 1,    3,
  "maize", list("Maize"),                        "Caod",  "cAOD",  0.01, 1
)

plots <- pmap(
  plot_specs,
  function(crop_name, crop_filter, var, x_label, hist_step, hist_scale) {
    
    adata <- data_all %>%
      filter(crop %in% unlist(crop_filter))
    
    p <- plot_interaction_response(
      adata = adata,
      var = var,
      x_label = x_label,
      subtitle_text = paste0(x_label, ": Faod × AOT40 vs no interaction"),
      hist_step = hist_step,
      hist_scale = hist_scale
    )
    
    ggsave(
      file.path(
        figure_dir,
        paste0("figureSI23_", tolower(var), "_interaction_comparison_", crop_name, ".tif")
      ),
      p,
      width = 6,
      height = 4.5,
      dpi = 300,
      device = "tiff",
      compression = "lzw"
    )
    
    p
  }
)

names(plots) <- paste(
  plot_specs$var,
  plot_specs$crop_name,
  sep = "_"
)
