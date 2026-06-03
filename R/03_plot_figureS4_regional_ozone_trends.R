source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

# Define TIFF folder and year range
tiff_folder <- "data/raw/ozone/O3_tif"
years <- 2001:2023

# Output directory
output_folder <- "figures/supplementary/figureS4"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# Load Chinese province boundaries and classify regions
province_plot <- st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json") %>%
  mutate(
    region = fcase(
      name %in% c("辽宁省", "吉林省", "黑龙江省"), "Northeast China",
      name %in% c("上海市", "江苏省", "浙江省", "安徽省", "福建省", "江西省", "山东省", "台湾省"), "East China",
      name %in% c("北京市", "天津市", "河北省", "山西省", "内蒙古自治区"), "North China",
      name %in% c("河南省", "湖北省", "湖南省"), "Central China",
      name %in% c("广东省", "广西壮族自治区", "海南省", "香港特别行政区", "澳门特别行政区"), "South China",
      name %in% c("四川省", "贵州省", "云南省", "西藏自治区", "重庆市"), "Southwest China",
      name %in% c("陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区"), "Northwest China"
    )
  )

# Aggregate provinces into regional polygons
region_plot <- province_plot %>%
  st_make_valid() %>%
  group_by(region) %>%
  summarise(geometry = st_union(geometry))

# Store regional annual means
region_data_list <- list()

# Extract regional mean O3 for each year
for (year in years) {
  
  tiff_path <- file.path(tiff_folder, paste0("south", year, ".tif"))
  
  if (!file.exists(tiff_path)) {
    warning(paste("File does not exist:", tiff_path))
    next
  }
  
  o3_raster <- rast(tiff_path)
  
  year_mean_values <- exact_extract(
    o3_raster,
    region_plot,
    "mean"
  )
  
  year_mean_df <- data.frame(
    region = region_plot$region,
    year = year,
    mean_O3 = year_mean_values
  )
  
  for (region_name in unique(year_mean_df$region)) {
    
    region_year_df <- year_mean_df %>%
      filter(region == region_name)
    
    if (!region_name %in% names(region_data_list)) {
      region_data_list[[region_name]] <- data.frame()
    }
    
    region_data_list[[region_name]] <- bind_rows(
      region_data_list[[region_name]],
      region_year_df
    )
  }
}

# Check one example region
print(region_data_list[["North China"]])

# Plot one example region
region_name <- "Northwest China"
region_data <- region_data_list[[region_name]]

breakpoint_year <- 2014

left_model <- lm(
  mean_O3 ~ year,
  data = region_data %>% filter(year <= breakpoint_year)
)

right_model <- lm(
  mean_O3 ~ year,
  data = region_data %>% filter(year > breakpoint_year)
)

left_end_y <- predict(
  left_model,
  newdata = data.frame(year = breakpoint_year)
)

right_start_y <- predict(
  right_model,
  newdata = data.frame(year = breakpoint_year)
)

p <- ggplot(region_data, aes(x = year, y = mean_O3)) +
  geom_line(color = "grey50", size = 0.8) +
  geom_segment(
    aes(
      x = min(region_data$year),
      y = coef(left_model)[1] + coef(left_model)[2] * min(region_data$year),
      xend = breakpoint_year,
      yend = left_end_y
    ),
    color = "red",
    size = 1
  ) +
  geom_segment(
    aes(
      x = breakpoint_year,
      y = right_start_y,
      xend = max(region_data$year),
      yend = coef(right_model)[1] + coef(right_model)[2] * max(region_data$year)
    ),
    color = "red",
    size = 1
  ) +
  geom_vline(
    xintercept = breakpoint_year,
    linetype = "dashed",
    color = "red",
    size = 0.9
  ) +
  scale_x_continuous(
    breaks = seq(min(region_data$year), max(region_data$year), by = 4),
    labels = seq(min(region_data$year), max(region_data$year), by = 4),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    limits = c(60, 110),
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = region_name,
    x = "year",
    y = "O3"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )

print(p)

# Remove empty or invalid regions
region_data_list <- region_data_list[
  !sapply(
    region_data_list,
    function(x) is.null(x) || nrow(x) == 0
  )
]

breakpoint_year <- 2014

# Export regional trend plots
for (region_name in names(region_data_list)) {
  
  region_data <- region_data_list[[region_name]]
  
  if (is.null(region_data) || nrow(region_data) == 0) {
    warning(paste("Region", region_name, "has no data and was skipped."))
    next
  }
  
  left_model <- lm(
    mean_O3 ~ year,
    data = region_data %>% filter(year <= breakpoint_year)
  )
  
  right_model <- lm(
    mean_O3 ~ year,
    data = region_data %>% filter(year > breakpoint_year)
  )
  
  left_end_y <- predict(
    left_model,
    newdata = data.frame(year = breakpoint_year)
  )
  
  right_start_y <- predict(
    right_model,
    newdata = data.frame(year = breakpoint_year)
  )
  
  p <- ggplot(region_data, aes(x = year, y = mean_O3)) +
    geom_line(color = "grey50", size = 0.8) +
    geom_segment(
      aes(
        x = min(region_data$year),
        y = coef(left_model)[1] + coef(left_model)[2] * min(region_data$year),
        xend = breakpoint_year,
        yend = left_end_y
      ),
      color = "red",
      size = 1
    ) +
    geom_segment(
      aes(
        x = breakpoint_year,
        y = right_start_y,
        xend = max(region_data$year),
        yend = coef(right_model)[1] + coef(right_model)[2] * max(region_data$year)
      ),
      color = "red",
      size = 1
    ) +
    geom_vline(
      xintercept = breakpoint_year,
      linetype = "dashed",
      color = "red",
      size = 0.9
    ) +
    scale_x_continuous(
      breaks = seq(min(region_data$year), max(region_data$year), by = 4),
      labels = seq(min(region_data$year), max(region_data$year), by = 4),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_y_continuous(
      limits = c(60, 110),
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    labs(
      title = region_name,
      x = "year",
      y = "O3"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 15, color = "black"),
      axis.text.y = element_text(size = 15, color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  tiff_filename <- file.path(
    output_folder,
    paste0(region_name, "_O3_Trend.tif")
  )
  
  tiff(
    tiff_filename,
    width = 2000,
    height = 1600,
    res = 300,
    compression = "lzw"
  )
  
  print(p)
  dev.off()
}