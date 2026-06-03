library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(readxl)

input_dir <- "data/derived/ablation"
figure_dir <- "figures/supplementary/ablation_interaction_gain"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Select one crop-specific ablation table
# df_ablation <- read_excel(file.path(input_dir, "ablation_maize.xlsx"))
# df_ablation <- read_excel(file.path(input_dir, "ablation_rice.xlsx"))
df_ablation <- read_excel(file.path(input_dir, "ablation_wheat.xlsx"))

factor_vars <- c(
  "Caod", "Faod", "AOT40", "cloud",
  "faodxAOT40",
  "bin", "root", "d", "co2"
)

plot_interaction_ablation_enhanced <- function(df, crop_name) {
  
  factor_vars <- c(
    "Caod", "Faod", "AOT40", "cloud",
    "faodxAOT40",
    "bin", "root", "d", "co2"
  )
  
  inter_var  <- "faodxAOT40"
  other_vars <- setdiff(factor_vars, inter_var)
  
  # Build matched model combinations
  df2 <- df %>%
    mutate(
      other_sum = rowSums(across(all_of(other_vars))),
      has_other = other_sum > 0
    ) %>%
    filter(has_other) %>%
    mutate(
      key = apply(across(all_of(other_vars)), 1, paste, collapse = "_"),
      label = apply(across(all_of(other_vars)), 1, function(z) {
        vars_on <- other_vars[z == 1]
        if (length(vars_on) == 0) "None" else paste(vars_on, collapse = " + ")
      })
    )
  
  # Keep combinations with both interaction and no-interaction versions
  df_pairs <- df2 %>%
    group_by(key) %>%
    filter(n_distinct(.data[[inter_var]]) == 2) %>%
    ungroup()
  
  # Average R² for each paired group
  df_pairs_mean <- df_pairs %>%
    mutate(
      interaction = ifelse(.data[[inter_var]] == 1, "with_int", "without_int")
    ) %>%
    group_by(key, label, interaction) %>%
    summarise(
      R2_mean = mean(.data[[crop_name]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Calculate R² gain from interaction
  df_diff <- df_pairs_mean %>%
    pivot_wider(
      names_from  = interaction,
      values_from = R2_mean
    ) %>%
    mutate(diff = with_int - without_int) %>%
    arrange(diff) %>%
    mutate(combo = factor(label, levels = label))
  
  df_sign <- df_diff %>%
    mutate(pos = diff > 0)
  
  # Prepare upset-style factor matrix
  df_factors <- df_pairs %>%
    filter(.data[[inter_var]] == 1) %>%
    group_by(key, label) %>%
    slice(1) %>%
    ungroup() %>%
    select(key, label, all_of(factor_vars)) %>%
    mutate(combo = factor(label, levels = levels(df_diff$combo)))
  
  df_long <- df_factors %>%
    select(combo, all_of(factor_vars)) %>%
    pivot_longer(
      cols      = all_of(factor_vars),
      names_to  = "Factor",
      values_to = "Included"
    )
  
  # Top panel: R² gain bar plot
  p_bar <- ggplot(df_diff, aes(x = combo, y = diff, fill = diff > 0)) +
    geom_col(color = "black") +
    scale_fill_manual(
      values = c(
        "TRUE" = "#F4A8A8",
        "FALSE" = "#9EC1E8"
      ),
      labels = c(
        "TRUE" = "R² increases",
        "FALSE" = "R² decreases"
      ),
      name = "Effect"
    ) +
    labs(
      y = paste0(crop_name, " R² gain (with − without interaction)"),
      x = NULL,
      title = crop_name
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "left",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Middle panel: sign indicator
  p_sign <- ggplot(df_sign, aes(x = combo, y = 1, color = pos)) +
    geom_point(size = 4) +
    scale_color_manual(
      values = c(
        "TRUE" = "#F4A8A8",
        "FALSE" = "#9EC1E8"
      ),
      guide = "none"
    ) +
    ylim(0.8, 1.2) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  # Bottom panel: upset-style inclusion matrix
  p_dot <- ggplot(df_long, aes(x = combo, y = Factor)) +
    geom_point(aes(color = Included == 1), size = 1.4) +
    scale_color_manual(
      values = c(
        "TRUE" = "black",
        "FALSE" = "#D3D3D3"
      ),
      guide = "none"
    ) +
    scale_y_discrete(limits = rev(factor_vars)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  p_bar / p_sign / p_dot +
    plot_layout(heights = c(3, 0.4, 2))
}

# Generate one crop figure
# p_maize <- plot_interaction_ablation_enhanced(df_ablation, "Maize")
# p_rice  <- plot_interaction_ablation_enhanced(df_ablation, "Rice")
p_wheat <- plot_interaction_ablation_enhanced(df_ablation, "Wheat")

# print(p_maize)
# print(p_rice)
print(p_wheat)

# Save one crop figure
# ggsave(file.path(figure_dir, "figureS13_interaction_ablation_gain_maize.tif"),
#        plot = p_maize,
#        width = 16, height = 8,
#        dpi = 600,
#        compression = "lzw")

# ggsave(file.path(figure_dir, "figureS14_interaction_ablation_gain_rice.tif"),
#        plot = p_rice,
#        width = 16, height = 8,
#        dpi = 600,
#        compression = "lzw")

ggsave(
  file.path(figure_dir, "figureS15_interaction_ablation_gain_wheat.tif"),
  plot = p_wheat,
  width = 16,
  height = 8,
  dpi = 600,
  compression = "lzw"
)