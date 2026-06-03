source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R") 

# Load wheat data
data <- get_data()
data <- data %>% filter(crop == "Wheat") %>% drop_na()

# Define models to fit
model_names <- c(
  "base" = "fml_base",
  "no_interaction" = "fml_base_no_interaction"
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
    fml_base, 
    fml_base_no_interaction
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
if(!dir.exists("figures/supplementary/si8")){
  dir.create("figures/supplementary/si8", recursive = TRUE)
}

ggsave("figures/supplementary/si8/response_climate_caod_Wheat.tif",
       plots,
       width=6,
       height=4.5,
       dpi=300,
       device="tiff")