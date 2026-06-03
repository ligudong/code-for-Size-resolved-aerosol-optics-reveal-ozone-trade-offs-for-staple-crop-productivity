source("R/00_load_packages.R")
source("R/00_helper_functions.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(raster)

crop_mask_dir <- "data/raw/crop_masks"
figure_dir <- "figures/figure1"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

rice <- raster(file.path(crop_mask_dir, "mask_Rice.tif"))
wheat <- raster(file.path(crop_mask_dir, "mask_Wheat.tif"))
maize <- raster(file.path(crop_mask_dir, "mask_Maize.tif"))

rice_df <- as.data.frame(rice, xy = TRUE)
wheat_df <- as.data.frame(wheat, xy = TRUE)
maize_df <- as.data.frame(maize, xy = TRUE)

names(rice_df)[3] <- "Rice"
names(wheat_df)[3] <- "Wheat"
names(maize_df)[3] <- "Maize"

crop_distribution <- rice_df %>%
  full_join(wheat_df, by = c("x", "y")) %>%
  full_join(maize_df, by = c("x", "y")) %>%
  pivot_longer(
    cols = c(Rice, Wheat, Maize),
    names_to = "crop",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

assign_region <- function(province_name) {
  fcase(
    province_name %in% c("辽宁省", "吉林省", "黑龙江省"), "Northeast China",
    province_name %in% c("上海市", "江苏省", "浙江省", "安徽省", "福建省", "江西省", "山东省", "台湾省"), "East China",
    province_name %in% c("北京市", "天津市", "河北省", "山西省", "内蒙古自治区"), "North China",
    province_name %in% c("河南省", "湖北省", "湖南省"), "Central China",
    province_name %in% c("广东省", "广西壮族自治区", "海南省", "香港特别行政区", "澳门特别行政区"), "South China",
    province_name %in% c("四川省", "贵州省", "云南省", "西藏自治区", "重庆市"), "Southwest China",
    province_name %in% c("陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区"), "Northwest China"
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

crop_plot <- ggplot() +
  geom_sf(
    data = region_plot,
    aes(fill = region),
    color = "black",
    alpha = 0.25,
    linewidth = 0.25,
    show.legend = FALSE
  ) +
  geom_tile(
    data = crop_distribution,
    aes(x = x, y = y, fill = crop)
  ) +
  scale_fill_manual(
    values = c(
      "Rice" = "#0050B350",
      "Maize" = "#FFA94050",
      "Wheat" = "#D4880680"
    ),
    na.value = "#D6E4FF"
  ) +
  labs(
    title = "Crop distribution in China (2019)",
    fill = "Crop"
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank()
  )

print(crop_plot)

ggsave(
  filename = file.path(figure_dir, "figure1c_crop_distribution_2019.tif"),
  plot = crop_plot,
  device = "tiff",
  dpi = 300,
  width = 8,
  height = 6,
  bg = "transparent",
  compression = "lzw"
)