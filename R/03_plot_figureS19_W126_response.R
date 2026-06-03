################fAOD###########
source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

figure_dir <- "figures/supplementary/si19"
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
    base = fml_base_W126,
    no_interaction = fml_W126_no_interaction
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
  file.path(figure_dir, "si19_faod_response_maize.tif"),
  plot_faod,
  width = 6,
  height = 4.5,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

ggsave(
  file.path(figure_dir, "si19_faod_response_maize.pdf"),
  plot_faod,
  width = 6,
  height = 4.5
)

##################cAOD###################
source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R") 

# Load wheat data
data <- get_data()
data <- data %>% filter(crop == "Wheat") %>% drop_na()

# Define models to fit
model_names <- c(
  "base" = "fml_base_W126",
  "no_interaction" = "fml_W126_no_interaction"
)

# Function to generate cAOD response plot
get_plot <- function(adata) {
  
  message("Processing crop data")
  message("Caod summary: ", summary(adata$Caod))
  
  data_modeling <- na.omit(adata)
  
  Caod_range <- range(data_modeling$Caod, na.rm = TRUE)
  if (!is.finite(Caod_range[1]) || !is.finite(Caod_range[2]) || Caod_range[1] == Caod_range[2]) {
    warning("Invalid Caod range, skip plot")
    return(NULL)
  }
  
  models <- lst(
    fml_base_W126, 
    fml_W126_no_interaction
  ) %>%
    map(~ feols(.x, data_modeling, cluster = ~city, weights = ~fraction, nthreads = 0, lean = TRUE), .progress = TRUE) %>%
    set_names(names(model_names))
  
  # Define x-axis for prediction
  Caod_hist <- seq(Caod_range[1], Caod_range[2], 0.001)
  Caod_mid  <- round(quantile(adata$Caod, 0.5, na.rm=TRUE),2)
  
  Caod_pred <- map_dfr(models, function(amodel) {
    relpred(amodel,
            tibble(Caod = Caod_hist),
            tibble(Caod = rep(Caod_mid, length(Caod_hist)))
    ) %>%
      mutate(Caod_x = Caod_hist)
  }, .id = "model") %>%
    mutate(across(c(fit,lwr,upr), ~ expm1(.x) * 1e2))
  
  Caod_offset <- (Caod_range[2]-Caod_range[1])*0.01
  
  # Labels for annotation
  Caod_lab <- map_dfr(models, function(amodel) {
    relpred(amodel,
            tibble(Caod = rep(Caod_range[2],10)),
            tibble(Caod = rep(Caod_mid, 10))
    ) %>%
      mutate(Caod_x = rep(Caod_range[2],10))
  }, .id = "model") %>%
    group_by(model) %>%
    slice_head(n=1) %>%
    ungroup() %>%
    mutate(across(c(fit,lwr,upr), ~ expm1(.x)*1e2)) %>%
    mutate(Caod_x = Caod_x + Caod_offset)
  
  # Define colors
  model_colors <- c(
    "base" = "#1f77b4",
    "no_interaction" = "#ff7f0e"
  )
  
  # Plot
  p1 <- Caod_pred %>%
    mutate(model = fct_relevel(model, names(model_colors))) %>%
    ggplot(aes(x = Caod_x, y = fit, color = model)) +
    geom_ribbon(aes(x = Caod_x, ymin=lwr, ymax=upr, fill=model),
                inherit.aes = FALSE,
                alpha = 0.2,
                data = . %>% filter(model=="base")) +
    geom_line(data = . %>% filter(model=="base"), size=1.5, linetype="solid") +
    geom_line(data = . %>% filter(model=="no_interaction"), size=1.5, linetype="solid") +
    scale_color_manual(values=model_colors) +
    geom_text_repel(aes(label=model),
                    data = Caod_lab,
                    show.legend = FALSE,
                    min.segment.length=0,
                    box.padding=1,
                    direction="y",
                    nudge_x=0.02,
                    nudge_y=0.02,
                    xlim=c(Caod_range[2]+Caod_offset, Caod_range[2]*1.1),
                    max.overlaps=100,
                    family="Roboto Condensed",
                    size=5,
                    seed=2022) +
    ylab("Percentage Change in SIF") +
    xlab("cAOD") +
    theme_half_open(18, font_family="Roboto Condensed") +
    background_grid() +
    theme(
      ggside.panel.scale=0.2,
      legend.title=element_blank(),
      legend.position=c(1.2,0.5),
      panel.grid=element_blank(),
      axis.text.x=element_text(size=50),
      axis.text.y=element_text(size=50)
    ) +
    geom_xsidehistogram(aes(x=Caod),
                        data=adata,
                        bins=80,
                        fill="#D3D3D3",
                        color="black",
                        size=0.2,
                        inherit.aes=FALSE) +
    ggside(x.pos="bottom", collapse="x") +
    scale_xsidey_continuous(labels=NULL, breaks=NULL) +
    scale_x_continuous(limits=c(Caod_range[1], Caod_range[2]+Caod_offset),
                       expand=expansion(mult=c(0,0.15))) +
    labs(subtitle="Crop Response to cAOD")
  
  return(p1)
}

# Generate plot
plots <- get_plot(data)

# Save
if(!dir.exists("figures/supplementary/si19")){
  dir.create("figures/supplementary/si19", recursive = TRUE)
}

ggsave("figures/supplementary/si19/response_climate_caod_Wheat.tif",
       plots,
       width=6,
       height=4.5,
       dpi=300,
       device="tiff")


###################W126#########################

source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")  # contains fml_base_W126 and fml_W126_no_interaction

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

# Function to generate W126 response plot
get_W126_response_plot <- function(adata) {
  
  message("Processing crop data")
  message("W126 summary: ", summary(adata$W126))
  
  data_modeling <- na.omit(adata)
  
  # Check range
  W126_range <- range(data_modeling$W126, na.rm = TRUE)
  if (!is.finite(W126_range[1]) || !is.finite(W126_range[2]) || W126_range[1]==W126_range[2]) {
    warning("Invalid W126 range, skip plot")
    return(NULL)
  }
  
  # Fit models
  models <- lst(
    fml_base_W126,
    fml_W126_no_interaction
  ) %>%
    map(~ feols(.x, data_modeling, cluster = ~city, weights = ~fraction, nthreads = 0, lean = TRUE),
        .progress = TRUE
    ) %>%
    set_names(model_names)
  
  # Prediction axis
  W126_hist <- seq(W126_range[1], W126_range[2], 1)
  W126_mid  <- round(quantile(adata$W126, 0.5, na.rm=TRUE), 2)
  W126_offset <- (W126_range[2] - W126_range[1]) * 0.01
  
  # Calculate predicted values
  W126_df <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(W126 = W126_hist),
      tibble(W126 = rep(W126_mid, length(W126_hist)))
    ) %>%
      mutate(W126_x = W126_hist)
  }, .id = "model") %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x)*1e2))
  
  # Labels for plotting
  W126_lab <- map_dfr(models, function(amodel) {
    relpred(
      amodel,
      tibble(W126 = rep(W126_range[2], 10)),
      tibble(W126 = rep(W126_mid, 10))
    ) %>%
      mutate(W126_x = rep(W126_range[2], 10))
  }, .id = "model") %>%
    group_by(model) %>%
    slice_head(n=1) %>%
    ungroup() %>%
    mutate(across(c(fit, lwr, upr), ~ expm1(.x)*1e2)) %>%
    mutate(W126_x = W126_x + W126_offset)
  
  # Plot
  model_colors <- c("base"="#1f77b4","no_interaction"="#ff7f0e")
  
  p1 <- W126_df %>%
    mutate(model = fct_relevel(model, rev(model_names))) %>%
    ggplot(aes(x=W126_x, y=fit, color=model)) +
    geom_ribbon(aes(x=W126_x, ymin=lwr, ymax=upr, fill=model),
                inherit.aes=FALSE, alpha=0.2, show.legend=FALSE,
                size=0.01, data=. %>% filter(model=="base")) +
    geom_line(data=. %>% filter(model=="base"), size=1.5) +
    geom_line(data=. %>% filter(model=="no_interaction"), size=1.5) +
    scale_color_manual(values=model_colors) +
    geom_text_repel(
      data=W126_lab,
      aes(label=model),
      min.segment.length=0,
      box.padding=1,
      direction="y",
      nudge_x=0.02,
      nudge_y=0.02,
      xlim=c(W126_range[2]+W126_offset, W126_range[2]*1.1),
      max.overlaps=100,
      family="Roboto Condensed",
      size=5,
      seed=2022
    ) +
    ylab("Percentage Change in SIF") +
    xlab("W126") +
    theme_half_open(18, font_family="Roboto Condensed") +
    background_grid() +
    theme(
      ggside.panel.scale=0.2,
      legend.title=element_blank(),
      legend.position=c(1.2,0.5),
      panel.grid=element_blank()
    ) +
    geom_xsidehistogram(aes(x=W126),
                        data=adata,
                        bins=80,
                        fill="#D3D3D3",
                        color="black",
                        size=0.2,
                        inherit.aes=FALSE) +
    ggside(x.pos="bottom", collapse="x") +
    scale_xsidey_continuous(labels=NULL, breaks=NULL) +
    scale_x_continuous(limits=c(W126_range[1], W126_range[2]+W126_offset),
                       expand=expansion(mult=c(0,0.15))) +
    labs(subtitle="Crop Response to W126")
  
  return(p1)
}

# Generate plots per crop
plot_maize <- get_W126_response_plot(data %>% filter(crop=="Maize"))
plot_rice  <- get_W126_response_plot(data %>% filter(crop %in% c("Rice(LR)","Rice(SR&ER)")))
plot_wheat <- get_W126_response_plot(data %>% filter(crop=="Wheat"))

# Save TIFF
ggsave(file.path(figure_dir,"si19_W126_response_maize.tif"), plot_maize, width=6, height=4.5, dpi=300, device="tiff")
ggsave(file.path(figure_dir,"si19_W126_response_rice.tif"),  plot_rice,  width=6, height=4.5, dpi=300, device="tiff")
ggsave(file.path(figure_dir,"si19_W126_response_wheat.tif"), plot_wheat, width=6, height=4.5, dpi=300, device="tiff")