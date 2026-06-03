library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(readxl)

input_dir <- "data/derived/ablation"
figure_dir <- "figures/supplementary/ablation_model_optimization"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

df_ablation <- read_excel(file.path(input_dir, "ablation_wheat.xlsx"))
#df_ablation <- read_excel(file.path(input_dir, "ablation_Maize.xlsx"))
#df_ablation <- read_excel(file.path(input_dir, "ablation_Rice.xlsx"))

# Factor variables for plotting (Upset plot)
factor_vars <- c("Caod", "Faod", "AOT40", "cloud", "faodxAOT40", "bin", "root", "d","co2")

# Set fml order for plotting
df_ablation <- df_ablation %>%
  mutate(fml_id = factor(fml_id, levels = paste0("fml_", 1:400)))

# Define plotting function
plot_crop_r2_ablation <- function(df, crop_name) {
  
  # Top panel: R² bar plot with colored bars for D and CO2 combinations
  p1 <- df %>%
    mutate(
      d_flag = d == 1,
      co2_flag = co2 == 1,
      dco2_group = case_when(
        d_flag & co2_flag ~ "both",
        d_flag & !co2_flag ~ "only_d",
        !d_flag & co2_flag ~ "only_co2",
        TRUE ~ "none"
      )
    ) %>%
    ggplot(aes(x = fml_id, y = .data[[crop_name]], fill = dco2_group)) +
    geom_col(width = 0.8, color = "black") +
    scale_fill_manual(
      values = c(
        "both" = "#F47F1E",
        "only_d" = "#6A8EAE90",
        "only_co2" = "#F3C80D",
        "none" = "white"
      ),
      name = "Model Extension"
    ) +
    labs(
      y = paste0(crop_name, " within R²"),
      x = NULL,
      title = crop_name
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Bottom panel: Upset-style scatter for factor inclusion
  df_long <- df %>%
    select(fml_id, all_of(factor_vars)) %>%
    pivot_longer(cols = all_of(factor_vars), names_to = "Factor", values_to = "Included") %>%
    mutate(color_group = ifelse(Included == 1, "included", "excluded"))
  
  p2 <- ggplot(df_long, aes(x = fml_id, y = Factor)) +
    geom_point(aes(color = color_group), size = 2) +
    scale_color_manual(values = c("included" = "black", "excluded" = "#D3D3D3")) +
    scale_y_discrete(limits = rev(factor_vars)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
      plot.margin = margin(0, 5, 5, 5),
      legend.position = "none"
    ) +
    labs(x = "Model (fml_1 ~ fml_400)", y = NULL)
  
  # Combine the bar plot and Upset plot vertically
  p <- p1 / p2 + plot_layout(heights = c(3, 1))
  
  return(p)
}

# Generate figure for Wheat (example)
p_wheat <- plot_crop_r2_ablation(df_ablation, "Wheat")
print(p_wheat)

ggsave(
  file.path(figure_dir, "figureS16_S18_interaction_ablation_wheat.tif"),
#  file.path(figure_dir, "figureS16_S18_interaction_ablation_Maize.tif"),
#  file.path(figure_dir, "figureS16_S18_interaction_ablation_Rice.tif"),
  plot = p_wheat,
  width = 30,
  height = 12,
  dpi = 600,
  compression = "lzw"
)