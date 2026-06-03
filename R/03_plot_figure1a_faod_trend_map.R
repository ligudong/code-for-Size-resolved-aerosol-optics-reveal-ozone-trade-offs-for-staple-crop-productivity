source("R/00_load_packages.R")
source("R/00_helper_functions.R")

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(exactextractr)
library(units)

# Directories
faod_dir <- "data/raw/faod_trend"
figure_dir <- "figures/figure1"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Function to assign regions based on province name
assign_region <- function(province_name) {
  fcase(
    province_name %in% c("辽宁省", "吉林省", "黑龙江省"), "Northeast China",
    province_name %in% c("上海市","江苏省","浙江省","安徽省","福建省","江西省","山东省","台湾省"), "East China",
    province_name %in% c("北京市","天津市","河北省","山西省","内蒙古自治区"), "North China",
    province_name %in% c("河南省","湖北省","湖南省"), "Central China",
    province_name %in% c("广东省","广西壮族自治区","海南省","香港特别行政区","澳门特别行政区"), "South China",
    province_name %in% c("四川省","贵州省","云南省","西藏自治区","重庆市"), "Southwest China",
    province_name %in% c("陕西省","甘肃省","青海省","宁夏回族自治区","新疆维吾尔自治区"), "Northwest China"
  )
}

# Load province boundaries
province_plot <- st_read(
  "https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json",
  quiet = TRUE
) %>%
  st_make_valid() %>%
  mutate(region = assign_region(name))

# Aggregate provinces into regions
region_plot <- province_plot %>%
  filter(!is.na(region)) %>%
  group_by(region) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# Load FAOD TIFF files by year
years <- 2001:2022
faod_files <- list.files(faod_dir, pattern = "\\.tif$", full.names = TRUE)
faod_stack <- rast(faod_files)
names(faod_stack) <- as.character(years)

# Calculate linear trend (slope) of FAOD for each pixel
faod_trend <- app(faod_stack, fun = function(x) {
  if(all(is.na(x))) return(NA_real_)
  coef(lm(x ~ years, na.action = na.exclude))[2]
})

# Project trend raster to province CRS
faod_trend <- project(faod_trend, crs(province_plot))

# Extract city-level mean trends
extracted <- exact_extract(faod_trend, province_plot, "mean") %>%
  bind_cols(province_plot, .) %>%
  rename(Faod_trend = mean) %>%
  mutate(area = units::drop_units(st_area(.))/1e6)

# Plot FAOD trend map
faod_trend_map <- ggplot(extracted) +
  geom_sf(aes(fill = Faod_trend), color = NA) +
  geom_sf(data = region_plot, fill = NA, color = "black", linewidth = 1.0) +
  scale_fill_gradientn(
    name = "fAOD Trend",
    colours = c("#20223e","#3b1f46","#7f4b89","#b46db3","#e3a5d6","#d8c2cb","#fcc893","#feb424","#fd8700"),
    values = scales::rescale(c(-0.02,-0.015,-0.01,-0.005,-0.002,0,0.005,0.01,0.015), to = c(0,1)),
    limits = c(-0.02,0.015),
    oob = scales::squish,
    na.value = "grey90",
    guide = guide_colorbar(title.position = "top", barwidth = 10, label.position = "bottom")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.2,0.2),
    legend.direction = "horizontal",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

print(faod_trend_map)

ggsave(
  file.path(figure_dir,"figure1a_faod_trend_map.tif"),
  plot = faod_trend_map,
  width = 10, height = 8, units = "in", dpi = 600,
  device = "tiff", compression = "lzw"
)

# -----------------------------
# Highlight pixels with p-value < 0.1
# -----------------------------

# Assume faod_p_value raster exists
faod_pvalue_file <- file.path(faod_dir, "faod_pvalue.tif")
if(file.exists(faod_pvalue_file)){
  faod_pvalue <- rast(faod_pvalue_file)
  faod_pvalue[faod_pvalue <= 0.1 | is.na(faod_pvalue)] <- NA
  faod_pvalue[faod_pvalue > 0.1] <- 1
  pvalue_points <- as.data.frame(faod_pvalue, xy = TRUE) %>%
    filter(!is.na(layer))
  
  pvalue_map <- ggplot() +
    geom_sf(data = region_plot, fill = NA, color = "black", linewidth = 1.0) +
    geom_point(data = pvalue_points, aes(x = x, y = y), color = "red", size = 0.8) +
    labs(title = "Pixels with FAOD trend p-value < 0.1") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  
  print(pvalue_map)
  
  ggsave(
    file.path(figure_dir,"figure1a_faod_pvalue_red_points.tif"),
    plot = pvalue_map,
    width = 10, height = 8, units = "in", dpi = 600,
    device = "tiff", compression = "lzw"
  )
}