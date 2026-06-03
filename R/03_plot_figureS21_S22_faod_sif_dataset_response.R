##########################fAOD###########################
source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

figure_dir <- "figures/supplementary/sif_robustness"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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

data <- get_data() %>%
  drop_na() %>%
  filter(crop %in% c("Rice(LR)", "Rice(SR&ER)"))
# filter(crop == "Maize")
# filter(crop == "Wheat")
plot_faod_response <- function(adata, formula_list, subtitle_text) {
  
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
  
  models <- formula_list %>%
    map(
      ~ feols(
        .x,
        data = data_modeling,
        cluster = ~ city,
        weights = ~ fraction,
        nthreads = 0,
        lean = FALSE
      )
    ) %>%
    set_names(model_names)
  
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
    ggplot(
      aes(
        x = Faod_x,
        y = fit,
        color = model,
        fill = model
      )
    ) +
    geom_ribbon(
      aes(ymin = lwr, ymax = upr),
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(
      linewidth = 1.5,
      show.legend = TRUE
    ) +
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
    ylab("Percentage change in SIF") +
    xlab("fAOD") +
    labs(subtitle = subtitle_text) +
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

plot_s21_faod <- plot_faod_response(
  data,
  formula_list = list(
    fml_base,
    fml_base2,
    fml_base3
  ),
  subtitle_text = "Response of GOSIF, CSIF and RTSIF to fAOD"
)

plot_s22_faod <- plot_faod_response(
  data,
  formula_list = list(
    fml_base_no_interaction,
    fml_base_no_interaction2,
    fml_base_no_interaction3
  ),
  subtitle_text = "Response of GOSIF, CSIF and RTSIF to fAOD without interaction"
)

print(plot_s21_faod)
print(plot_s22_faod)

ggsave(
  file.path(figure_dir, "figureS21_faod_sif_dataset_response_rice.tif"),
  plot = plot_s21_faod,
  width = 6.5,
  height = 4.8,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(figure_dir, "figureS22_faod_sif_dataset_response_rice.tif"),
  plot = plot_s22_faod,
  width = 6.5,
  height = 4.8,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)


#####################cAOD########################
source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")  # fml_base, fml_base2, fml_base3 / fml_base_no_interaction, etc.


# Output directory
figure_dir <- "figures/supplementary/sif_robustness"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
data <- get_data() %>% drop_na()

# Filter crop (Rice combined)
data <- data %>% filter(crop %in% c("Rice(LR)", "Rice(SR&ER)"))

# Model names
model_names <- c("GOSIF","CSIF","RTSIF")

get_plot <- function(adata, formula_list, subtitle_text) {
  
  message("Processing crop data")
  message("Caod summary: ", summary(adata$Caod))
  
  data_modeling <- na.omit(adata)
  Caod_range <- range(data_modeling$Caod, na.rm = TRUE)
  if (!is.finite(Caod_range[1]) || !is.finite(Caod_range[2]) || Caod_range[1]==Caod_range[2]) {
    warning("Invalid Caod range, skip plot")
    return(NULL)
  }
  
  # Fit three models: GOSIF, CSIF, RTSIF
  models <- formula_list %>%
    map(~ feols(.x, data_modeling, cluster=~city, weights=~fraction, nthreads=0, lean=FALSE), .progress=TRUE) %>%
    set_names(model_names)
  
  # Construct Caod response axis
  Caod_hist <- seq(Caod_range[1], Caod_range[2], by=0.001)
  Caod_mid <- round(quantile(data_modeling$Caod,0.5,na.rm=TRUE),2)
  
  # Predicted values
  Caod_pred <- map_dfr(models, function(amodel) {
    relpred(amodel, tibble(Caod=Caod_hist),
            tibble(Caod=rep(Caod_mid, length(Caod_hist)))) %>%
      mutate(Caod_x=Caod_hist)
  }, .id="model") %>%
    mutate(across(c(fit,lwr,upr), ~ expm1(.x)*100))
  
  Caod_offset <- (Caod_range[2]-Caod_range[1])*0.02
  
  # Labels for curves
  Caod_lab <- map_dfr(models, function(amodel) {
    relpred(amodel, tibble(Caod=rep(Caod_range[2],10)),
            tibble(Caod=rep(Caod_mid,10))) %>%
      mutate(Caod_x=rep(Caod_range[2],10))
  }, .id="model") %>%
    group_by(model) %>% slice_head(n=1) %>% ungroup() %>%
    mutate(across(c(fit,lwr,upr), ~ expm1(.x)*100),
           Caod_x=Caod_x+Caod_offset)
  
  # Colors
  model_colors <- c(
    "GOSIF"="#1f77b4",
    "CSIF"=scales::alpha("#2ca02c",0.35),
    "RTSIF"=scales::alpha("#d62728",0.35)
  )
  
  model_fills <- c(
    "GOSIF"=scales::alpha("#1f77b4",0.15),
    "CSIF"=scales::alpha("#2ca02c",0.06),
    "RTSIF"=scales::alpha("#d62728",0.06)
  )
  
  p1 <- Caod_pred %>%
    mutate(model=fct_relevel(model,model_names)) %>%
    ggplot(aes(x=Caod_x, y=fit, color=model, fill=model)) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), color=NA, show.legend=FALSE) +
    geom_line(size=1.5) +
    geom_text_repel(aes(label=model), data=Caod_lab,
                    show.legend=FALSE, min.segment.length=0, box.padding=0.8,
                    direction="y", nudge_x=0.02, nudge_y=0.02,
                    max.overlaps=100, family="Roboto Condensed", size=5, seed=2022) +
    geom_xsidehistogram(aes(x=Caod), data=adata, bins=80,
                        fill="#D3D3D3", color="black", size=0.2, inherit.aes=FALSE) +
    scale_color_manual(values=model_colors) +
    scale_fill_manual(values=model_fills) +
    ylab("Percentage Change in SIF") +
    xlab("cAOD") +
    labs(subtitle=subtitle_text) +
    theme_half_open(18,font_family="Roboto Condensed") +
    background_grid() +
    ggside(x.pos="bottom", collapse="x") +
    scale_xsidey_continuous(labels=NULL, breaks=NULL) +
    scale_x_continuous(limits=c(Caod_range[1],Caod_range[2]+Caod_offset),
                       expand=expansion(mult=c(0,0.15))) +
    theme(ggside.panel.scale=0.45, legend.title=element_blank(),
          legend.position="right", panel.grid=element_blank())
  
  return(p1)
}

# S21: base interaction
plot_s21 <- get_plot(
  data,
  formula_list=list(fml_base,fml_base2,fml_base3),
  subtitle_text="Response of GOSIF, CSIF, RTSIF to cAOD (interaction)"
)

# S22: no interaction
plot_s22 <- get_plot(
  data,
  formula_list=list(fml_base_no_interaction,fml_base_no_interaction2,fml_base_no_interaction3),
  subtitle_text="Response of GOSIF, CSIF, RTSIF to cAOD (no interaction)"
)

# Print plots
print(plot_s21)
print(plot_s22)

# Save plots
ggsave(file.path(figure_dir,"figureS21_caod_sif_dataset_response_rice.tif"), plot_s21,
       width=6.5, height=4.8, dpi=300, device="tiff", compression="lzw")

ggsave(file.path(figure_dir,"figureS22_caod_sif_dataset_response_rice.tif"), plot_s22,
       width=6.5, height=4.8, dpi=300, device="tiff", compression="lzw")

###############################AOT40################################
source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

figure_dir <- "figures/supplementary/sif_robustness"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

data <- get_data() %>%
  drop_na()

data <- data %>%
  filter(crop == "Wheat")

model_names <- c("GOSIF", "CSIF", "RTSIF")

get_plot <- function(adata, formula_list, subtitle_text) {
  
  message("Processing crop data")
  message("AOT40 summary: ", summary(adata$AOT40))
  
  data_modeling <- na.omit(adata)
  
  AOT40_range <- range(data_modeling$AOT40, na.rm = TRUE)
  
  if (
    !is.finite(AOT40_range[1]) ||
    !is.finite(AOT40_range[2]) ||
    AOT40_range[1] == AOT40_range[2]
  ) {
    warning("Invalid AOT40 range, skip plot")
    return(NULL)
  }
  
  models <- formula_list %>%
    map(
      ~ feols(
        .x,
        data = data_modeling,
        cluster = ~ city,
        weights = ~ fraction,
        nthreads = 0,
        lean = FALSE
      ),
      .progress = TRUE
    ) %>%
    set_names(model_names)
  
  AOT40_hist <- seq(
    AOT40_range[1],
    AOT40_range[2],
    by = 1
  )
  
  AOT40_mid <- round(
    quantile(
      data_modeling$AOT40,
      0.5,
      na.rm = TRUE
    ),
    2
  )
  
  AOT40_pred <- map_dfr(
    models,
    function(amodel) {
      relpred(
        amodel,
        tibble(AOT40 = AOT40_hist),
        tibble(AOT40 = rep(AOT40_mid, length(AOT40_hist)))
      ) %>%
        mutate(AOT40_x = AOT40_hist)
    },
    .id = "model"
  ) %>%
    mutate(
      across(
        c(fit, lwr, upr),
        ~ expm1(.x) * 100
      )
    )
  
  AOT40_offset <- (AOT40_range[2] - AOT40_range[1]) * 0.02
  
  AOT40_lab <- map_dfr(
    models,
    function(amodel) {
      relpred(
        amodel,
        tibble(AOT40 = rep(AOT40_range[2], 10)),
        tibble(AOT40 = rep(AOT40_mid, 10))
      ) %>%
        mutate(AOT40_x = rep(AOT40_range[2], 10))
    },
    .id = "model"
  ) %>%
    group_by(model) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(
      across(
        c(fit, lwr, upr),
        ~ expm1(.x) * 100
      ),
      AOT40_x = AOT40_x + AOT40_offset
    )
  
  model_colors <- c(
    "GOSIF" = "#1f77b4",
    "CSIF" = scales::alpha("#2ca02c", 0.35),
    "RTSIF" = scales::alpha("#d62728", 0.35)
  )
  
  model_fills <- c(
    "GOSIF" = scales::alpha("#1f77b4", 0.15),
    "CSIF" = scales::alpha("#2ca02c", 0.06),
    "RTSIF" = scales::alpha("#d62728", 0.06)
  )
  
  p1 <- AOT40_pred %>%
    mutate(
      model = fct_relevel(
        model,
        model_names
      )
    ) %>%
    ggplot(
      aes(
        x = AOT40_x,
        y = fit,
        color = model,
        fill = model
      )
    ) +
    geom_ribbon(
      aes(
        ymin = lwr,
        ymax = upr
      ),
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(
      size = 1.5,
      show.legend = TRUE
    ) +
    geom_text_repel(
      aes(label = model),
      data = AOT40_lab,
      show.legend = FALSE,
      min.segment.length = 0,
      box.padding = 0.8,
      direction = "y",
      nudge_x = 0.02,
      nudge_y = 0.02,
      max.overlaps = 100,
      family = "Roboto Condensed",
      size = 5,
      seed = 2022
    ) +
    geom_xsidehistogram(
      aes(
        x = AOT40,
        y = after_stat(count * 3)
      ),
      data = adata,
      bins = 80,
      fill = "#D3D3D3",
      color = "black",
      size = 0.2,
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = model_colors
    ) +
    scale_fill_manual(
      values = model_fills
    ) +
    ylab("Percentage Change in SIF") +
    xlab("AOT40") +
    labs(
      subtitle = subtitle_text
    ) +
    theme_half_open(
      18,
      font_family = "Roboto Condensed"
    ) +
    background_grid() +
    ggside(
      x.pos = "bottom",
      collapse = "x"
    ) +
    scale_xsidey_continuous(
      labels = NULL,
      breaks = NULL
    ) +
    scale_x_continuous(
      limits = c(
        AOT40_range[1],
        AOT40_range[2] + AOT40_offset
      ),
      expand = expansion(
        mult = c(0, 0.15)
      )
    ) +
    theme(
      ggside.panel.scale = 3,
      legend.title = element_blank(),
      legend.position = "right",
      panel.grid = element_blank()
    )
  
  return(p1)
}

plot_s21 <- get_plot(
  data,
  formula_list = list(
    fml_base,
    fml_base2,
    fml_base3
  ),
  subtitle_text = "Response of GOSIF, CSIF, RTSIF to AOT40 (interaction)"
)

plot_s22 <- get_plot(
  data,
  formula_list = list(
    fml_base_no_interaction,
    fml_base_no_interaction2,
    fml_base_no_interaction3
  ),
  subtitle_text = "Response of GOSIF, CSIF, RTSIF to AOT40 (no interaction)"
)

print(plot_s21)
print(plot_s22)

ggsave(
  file.path(
    figure_dir,
    "figureS21_aot40_sif_dataset_response_wheat.tif"
  ),
  plot_s21,
  width = 6.5,
  height = 4.8,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(
    figure_dir,
    "figureS22_aot40_sif_dataset_response_wheat.tif"
  ),
  plot_s22,
  width = 6.5,
  height = 4.8,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)