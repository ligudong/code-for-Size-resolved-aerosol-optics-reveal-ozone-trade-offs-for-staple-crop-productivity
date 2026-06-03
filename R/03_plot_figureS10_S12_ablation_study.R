# Load necessary packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(readxl)

# Define input path for ablation study results
# Replace previous hard-coded paths with a structured folder
input_file <- "data/derived/ablation/ablation_results.xlsx"  

# Read the Excel file (replace crop-specific comment with unified path)
# Example: choose which sheet or file corresponds to Rice
df_ablation <- read_excel(input_file)

# Factor variables for plotting Upset plot
factor_vars <- c("Caod", "Faod", "AOT40", "cloud", "faodxAOT40", "bin", "root", "d", "co2")

# Set the formula order for plotting
df_ablation <- df_ablation %>%
  mutate(fml_id = factor(fml_id, levels = paste0("fml_", 1:400)))

# Define plotting function
plot_crop_r2_ablation <- function(df, crop_name) {
  # Extract R² values for this crop
  crop_r2 <- df %>% select(fml_id, all_of(crop_name))
  
  # Interaction column
  interact_var <- "faodxAOT40"
  
  # Top panel: R² bar chart
  p1 <- df %>% 
    mutate(is_interaction = ifelse(.data[[interact_var]] == 1, "interaction", "no_interaction")) %>%
    ggplot(aes(x = fml_id, 
               y = .data[[crop_name]],
               fill = is_interaction)) +
    geom_col(width = 0.8, color = "black") +
    scale_fill_manual(values = c("interaction" = "#FC9272", 
                                 "no_interaction" = "#6BAED660")) +
    labs(y = paste0(crop_name, " within R²"), x = NULL, fill = "Model Type") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      plot.margin = margin(5, 5, 0, 5)
    )
  
  # Bottom panel: Upset scatter plot
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
  
  # Combine panels
  p <- p1 / p2 + plot_layout(heights = c(3, 1))
  return(p)
}

# Generate the plot for Rice
p_rice <- plot_crop_r2_ablation(df_ablation, "Rice")

# Print plot
print(p_rice)

# Define output directory
output_dir <- "figures/supplementary/si10_si12"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save as high-resolution TIFF
ggsave(file.path(output_dir, "figureSI10_Rice_r2_ablation.tif"), 
       plot = p_rice, 
       width = 30, height = 12, 
       dpi = 600, 
       compression = "lzw")