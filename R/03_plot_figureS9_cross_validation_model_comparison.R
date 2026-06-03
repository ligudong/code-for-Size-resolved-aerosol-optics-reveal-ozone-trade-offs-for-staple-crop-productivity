source("R/00_load_packages.R")      
source("R/00_helper_functions.R")    
source("R/00_formula_definitions.R") 


data <- get_data()


n_folds <- 10

data <- data %>% nest(fdata = -crop_parent)

set.seed(i) 

formulas <- list(
  fml_1 = fml_faod,
  fml_2 = fml_aot40,
  fml_3 = fml_af,
  fml_4 = fml_faod_aot40_caod
)


data <- data %>%
  mutate(coefs = map(fdata, function(adata) {

    adata <- adata %>%
      nest(.by = city, .key = "ffdata")

    results <- map_dfr(1:n_folds, function(fold) {

      set.seed(fold)

      aadata <- adata %>%
        unnest(ffdata) %>%
        group_by(city, year) %>%
        mutate(across(c(Caod, Faod, cloud, AOT40, bin1:bin42, root_1:root_9, d1:d13), ~ . - mean(.), .names = "demeaned_{.col}")) %>%
        ungroup()
      
      map_dfr(formulas, function(fml) {
        
        model <- feols(fml, aadata,
                       weights = ~fraction,
                       nthreads = 0, notes = FALSE, lean = TRUE,
                       cluster = ~x_y)  # 按年份（x_y）聚类标准误差
        
        tidy_res <- tidy(model)
        conf_int <- confint(model)
        tidy_res <- cbind(tidy_res, conf_int)  # 合并置信区间到回归结果中
        
        glance_res <- glance(model)
        
        tidy_res %>%
          mutate(
            r_squared = glance_res$r.squared,  # R²
            adj_r_squared = glance_res$adj.r.squared,  # Adjusted R²
            within_r_squared = glance_res$within.r.squared,  # within R²
            fold_id = fold  # 标记当前折数
          )
        
      }, .id = "formula", .progress = TRUE)
      
    }, .id = "id", .progress = TRUE)
    
    results
  }))

results_summary <- data %>%
  mutate(
    coefs_summary = map(coefs, function(coef_data) {
      coef_data %>%
        group_by(term, formula) %>%
        summarise(
          mean_estimate = mean(estimate, na.rm = TRUE), 
          mean_p_value = mean(p.value, na.rm = TRUE), 
          mean_r_squared = mean(r_squared, na.rm = TRUE), 
          mean_adj_r_squared = mean(adj_r_squared, na.rm = TRUE), 
          mean_within_r_squared = mean(within_r_squared, na.rm = TRUE) 
        )
    })
  ) %>%
  select(crop_parent, coefs_summary)  

final_results <- results_summary %>%
  unnest(coefs_summary)


print(final_results)


final_results_filtered <- final_results %>%
  select(crop_parent, formula, mean_within_r_squared)


final_results_filtered <- final_results %>%
  select(crop_parent, formula, mean_within_r_squared) %>%
  distinct()


ggplot(final_results_filtered, aes(x = formula, y = mean_within_r_squared, color = crop_parent, shape = crop_parent, group = crop_parent)) +
  geom_point(size = 4, alpha = 0.7) + 
  geom_line(aes(linetype = crop_parent), size = 1) +  
  labs(title = "Mean within R-squared by Formula and Crop",
       x = "Formula",
       y = "Mean within R-squared") +
  theme_minimal() +
  scale_shape_manual(values = c(16, 17, 18)) +  # 使用不同的形状
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", size = 14),  
        axis.text.y = element_text(family = "Arial", size = 14),  
        axis.title.x = element_text(family = "Arial", size = 16),  
        axis.title.y = element_text(family = "Arial", size = 16), 
        plot.title = element_text(family = "Arial", size = 18, face = "bold"), 
        legend.title = element_text(family = "Arial", size = 14), 
        legend.text = element_text(family = "Arial", size = 12)) 



plot <- ggplot(final_results_filtered, aes(x = formula, y = mean_within_r_squared, color = crop_parent, shape = crop_parent, group = crop_parent)) +
  geom_point(size = 4, alpha = 0.7) +  
  geom_line(aes(linetype = crop_parent), size = 1) + 
  labs(
       x = "Formula",
       y = "Mean within R-squared") +
  theme_minimal() +
  scale_shape_manual(values = c(16, 17, 18)) +  # 使用不同的形状
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", size = 14), 
        axis.text.y = element_text(family = "Arial", size = 45),  
        axis.title.x = element_text(family = "Arial", size = 45),  
        axis.title.y = element_text(family = "Arial", size = 45),  
        plot.title = element_text(family = "Arial", size = 18, face = "bold"),  
        legend.title = element_text(family = "Arial", size = 40), 
        legend.text = element_text(family = "Arial", size = 40)) 


ggsave(filename = "figures/supplementary/figureS9/mean_within_r_squared_plot.tif", plot = plot, device = "tiff", dpi = 300, width = 10, height = 8)

write.csv(final_results_filtered,
          file = "figures/supplementary/figureS9/results_filtered.csv",
          row.names = FALSE)