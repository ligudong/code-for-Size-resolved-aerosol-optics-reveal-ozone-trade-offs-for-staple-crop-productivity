source("R/00_load_packages.R")
source("R/00_helper_functions.R")

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(exactextractr)
library(units)

ozone_dir <- "data/derived/ozone_yearly_tifs"
figure_dir <- "figures/figure1"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

assign_region <- function(province_name) {
  fcase(
    province_name %in% c("辽宁省","吉林省","黑龙江省"), "Northeast China",
    province_name %in% c("上海市","江苏省","浙江省","安徽省","福建省","江西省","山东省","台湾省"), "East China",
    province_name %in% c("北京市","天津市","河北省","山西省","内蒙古自治区"), "North China",
    province_name %in% c("河南省","湖北省","湖南省"), "Central China",
    province_name %in% c("广东省","广西壮族自治区","海南省","香港特别行政区","澳门特别行政区"), "South China",
    province_name %in% c("四川省","贵州省","云南省","西藏自治区","重庆市"), "Southwest China",
    province_name %in% c("陕西省","甘肃省","青海省","宁夏回族自治区","新疆维吾尔自治区"), "Northwest China"
  )
}

province_plot <- st_read(
  "https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json",
  quiet = TRUE
) %>%
  st_make_valid() %>%
  mutate(region = assign_region(name))

region_plot <- province_plot %>%
  filter(!is.na(region)) %>%
  group_by(region) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

years <- 2005:2022
ozone_files <- list.files(ozone_dir, pattern="O3_X\\d{4}\\.tif$", full.names = TRUE)
ozone_files <- ozone_files[order(readr::parse_number(basename(ozone_files)))]
ozone_stack <- rast(ozone_files)
names(ozone_stack) <- as.character(years)

# Compute trend and p-value
ozone_trend_with_p <- app(ozone_stack, fun=function(x){
  if(all(is.na(x)) || sum(!is.na(x)) < 3) return(c(NA_real_, NA_real_))
  model <- tryCatch(lm(x ~ years, na.action=na.exclude), error=function(e) NULL)
  if(is.null(model)) return(c(NA_real_, NA_real_))
  c(slope=coef(model)[2], p_value=summary(model)$coefficients[2,4])
})

ozone_trend <- project(ozone_trend_with_p[[1]], crs(province_plot))
ozone_pvalue <- project(ozone_trend_with_p[[2]], crs(province_plot))

# Extract city-level mean
extracted <- exact_extract(ozone_trend, province_plot, "mean") %>%
  bind_cols(province_plot, .) %>%
  rename(O3_trend = mean) %>%
  mutate(area = units::drop_units(st_area(.))/1e6)

# Plot trend map
trend_map <- ggplot(extracted) +
  geom_sf(aes(fill=O3_trend), color=NA) +
  geom_sf(data=region_plot, fill=NA, color="black", linewidth=1) +
  scale_fill_gradientn(
    name = expression(O[3]~trend~"(ppb yr"^{-1}*")"),
    colours = c("#585885","#8080ad","#9c9ccc","#cacae0","#f1d7da","#f0a8a6","#e8817d","#be5958"),
    limits=c(-2,2),
    oob=scales::squish,
    na.value="grey90"
  ) +
  theme_minimal(base_size=14) +
  theme(
    legend.position=c(0.2,0.2),
    legend.direction="horizontal",
    axis.text=element_blank(),
    axis.title=element_blank(),
    panel.grid=element_blank()
  )

ggsave(
  file.path(figure_dir,"figure1b_ozone_trend_map.tif"),
  plot=trend_map,
  width=10,
  height=8,
  units="in",
  dpi=600,
  device="tiff",
  compression="lzw"
)

# Mark pixels with p-value > 0.1 using black ×
significant_points <- as.data.frame(ozone_pvalue, xy=TRUE, na.rm=FALSE) %>%
  rename(O3_p_value = 3) %>%
  filter(!is.na(O3_p_value), O3_p_value > 0.1)

pvalue_map <- ggplot() +
  geom_sf(data=region_plot, fill=NA, color="black", linewidth=1) +
  geom_point(data=significant_points, aes(x=x, y=y), color="black", size=1.2, shape=4) +
  coord_sf(expand=FALSE) +
  labs(title=expression(O[3]~trend~italic(p)~"> 0.1")) +
  theme_minimal(base_size=14) +
  theme(
    axis.text=element_blank(),
    axis.title=element_blank(),
    panel.grid=element_blank(),
    plot.title=element_text(hjust=0.5, face="bold")
  )

ggsave(
  file.path(figure_dir,"figure1b_ozone_pvalue_black_x.tif"),
  plot=pvalue_map,
  width=10,
  height=8,
  units="in",
  dpi=600,
  device="tiff",
  compression="lzw"
)