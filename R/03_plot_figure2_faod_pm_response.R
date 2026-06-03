source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

# Directories
figure_dir <- "figures/figure2"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Load data for modeling (already cleaned)
data <- get_data() %>% drop_na()

# Filter crop if needed
data <- data %>% filter(crop == "Rice")  # change to "Maize" or "Wheat" as needed
data_modeling <- na.omit(data)

# Fit three models
models <- list(
  base = feols(
    fml_base,
    data = data_modeling,
    cluster = ~city,
    weights = ~fraction,
    nthreads = 0,
    lean = TRUE
  ),
  no_jiaohu = feols(
    fml_base_no_interaction,
    data = data_modeling,
    cluster = ~city,
    weights = ~fraction,
    nthreads = 0,
    lean = TRUE
  ),
  base_pm = feols(
    fml_pm_interaction_gosif,
    data = data_modeling,
    cluster = ~city,
    weights = ~fraction,
    nthreads = 0,
    lean = TRUE
  )
)

# Map PM <-> fAOD
faod_min <- 0
faod_max <- 1.5
pm_min <- 0
pm_max <- 80
b_map <- (faod_max - faod_min)/(pm_max - pm_min)
a_map <- faod_min - b_map * pm_min

pm_to_faod <- function(pm) a_map + b_map * pm
faod_to_pm <- function(faod) (faod - a_map)/b_map

# Construct fAOD curves
Faod_hist <- seq(faod_min, faod_max, length.out = 300)
Faod_mid <- round(quantile(data_modeling$Faod, 0.5, na.rm=TRUE), 2)

Faod_curve <- map_dfr(
  models[c("base", "no_jiaohu")],
  function(amodel) {
    relpred(amodel,
            tibble(Faod = Faod_hist),
            tibble(Faod = rep(Faod_mid, length(Faod_hist)))
    ) %>%
      mutate(Faod_x = Faod_hist)
  },
  .id = "model"
) %>%
  mutate(
    x_plot = Faod_x,
    across(c(fit,lwr,upr), ~ expm1(.x)*1e2)
  )

# Construct PM2.5 curve
PM_hist <- seq(pm_min, pm_max, length.out = 300)
PM_mid <- round(quantile(data_modeling$PM25, 0.5, na.rm=TRUE), 2)

PM_curve <- relpred(
  models$base_pm,
  tibble(PM25 = PM_hist),
  tibble(PM25 = rep(PM_mid, length(PM_hist)))
) %>%
  mutate(
    model = "base_pm",
    PM25_x = PM_hist,
    x_plot = pm_to_faod(PM25_x)
  ) %>%
  mutate(across(c(fit,lwr,upr), ~ expm1(.x)*1e2))

# Merge curves
curve_data <- bind_rows(
  Faod_curve %>% mutate(line_group = case_when(model=="base" ~ "fAOD base",
                                               model=="no_jiaohu" ~ "fAOD no interaction")),
  PM_curve %>% mutate(line_group="PM2.5 base")
)

# Histogram data
pm_hist_data <- data_modeling %>% filter(PM25 >= pm_min, PM25 <= pm_max) %>% transmute(x_plot=pm_to_faod(PM25))
faod_hist_data <- data_modeling %>% filter(Faod >= faod_min, Faod <= faod_max) %>% transmute(x_plot=Faod)

# Axis breaks
faod_breaks <- seq(0,1.5,0.3)
pm_breaks <- seq(0,80,20)
pm_break_positions <- pm_to_faod(pm_breaks)
x_limits <- c(0,1.5)

# Colors
line_colors <- c("fAOD base"="#2ca02c","fAOD no interaction"="#ff7f0e","PM2.5 base"="#1f77b4")
pm_hist_fill <- "#B8C4D6"
faod_hist_fill <- "#D8C3A5"

# Top: PM2.5 inverted histogram
p_top <- ggplot(pm_hist_data, aes(x=x_plot)) +
  geom_histogram(aes(y=-after_stat(count)), bins=55, fill=pm_hist_fill, color="black", linewidth=0.25) +
  scale_x_continuous(position="top", breaks=pm_break_positions, labels=pm_breaks, limits=x_limits, expand=c(0,0)) +
  labs(x=expression(PM[2.5]), y=NULL) +
  theme_classic(base_size=15)

# Middle: response curves
p_mid <- ggplot() +
  geom_hline(yintercept=0, color="grey35", linewidth=0.5, linetype="22") +
  geom_ribbon(data=curve_data %>% filter(line_group=="PM2.5 base"),
              aes(x=x_plot, ymin=lwr, ymax=upr), fill=line_colors["PM2.5 base"], alpha=0.12) +
  geom_ribbon(data=curve_data %>% filter(line_group=="fAOD base"),
              aes(x=x_plot, ymin=lwr, ymax=upr), fill=line_colors["fAOD base"], alpha=0.10) +
  geom_line(data=curve_data %>% filter(line_group=="PM2.5 base"), aes(x=x_plot, y=fit, color=line_group), linewidth=1.5) +
  geom_line(data=curve_data %>% filter(line_group=="fAOD no interaction"), aes(x=x_plot, y=fit, color=line_group), linewidth=1.4, linetype="42") +
  geom_line(data=curve_data %>% filter(line_group=="fAOD base"), aes(x=x_plot, y=fit, color=line_group), linewidth=1.5) +
  scale_color_manual(values=line_colors) +
  scale_x_continuous(limits=x_limits) +
  ylab("Percentage Change in SIF") +
  xlab(NULL) +
  theme_minimal(base_size=15) +
  theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank())

# Bottom: fAOD histogram
p_bottom <- ggplot(faod_hist_data, aes(x=x_plot)) +
  geom_histogram(bins=55, fill=faod_hist_fill, color="black", linewidth=0.25) +
  scale_x_continuous(breaks=faod_breaks, labels=faod_breaks, limits=x_limits, expand=c(0,0)) +
  labs(x="fAOD", y=NULL) +
  theme_classic(base_size=15)

# Combine plots
plot_final <- p_top / p_mid / p_bottom + plot_layout(heights=c(1,3.9,1))
print(plot_final)

# Save TIFF
dir.create("figures/crop/p2/2026", recursive=TRUE, showWarnings=FALSE)
tiff(filename="figures/crop/p2/2026/faod_pm_transparent_Rice.tif",
     width=8,height=9,units="in",res=300,compression="lzw",bg="transparent",type="cairo")
print(plot_final)
dev.off()

tiff(filename="figures/crop/p2/2026/faod_pm_opaque_Rice.tif",
     width=8,height=9,units="in",res=300,compression="lzw",bg="white",type="cairo")
print(plot_final)
dev.off()