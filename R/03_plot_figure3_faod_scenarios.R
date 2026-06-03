source("R/00_load_packages.R")
source("R/00_helper_functions.R")

library(sf)
library(terra)
library(qs)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(exactextractr)
library(writexl)

province_shp_file <- "data/raw/shp/province.shp"
impact_file <- "data/derived/impacts/impacts_Faod.qs"
summary_file <- "data/derived/impacts/impacts_Faod_summarised.rds"

figure_dir <- "figures/faod"
table_dir <- "data/derived/impacts"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

region_order <- c(
  "North China",
  "Northeast China",
  "East China",
  "Central China",
  "South China",
  "Southwest China",
  "Northwest China",
  "China"
)

assign_region <- function(province_name) {
  fcase(
    province_name %in% c("辽宁省", "吉林省", "黑龙江省"),
    "Northeast China",
    province_name %in% c("上海市", "江苏省", "浙江省", "安徽省", "福建省", "江西省", "山东省", "台湾省"),
    "East China",
    province_name %in% c("北京市", "天津市", "河北省", "山西省", "内蒙古自治区"),
    "North China",
    province_name %in% c("河南省", "湖北省", "湖南省"),
    "Central China",
    province_name %in% c("广东省", "广西壮族自治区", "海南省", "香港特别行政区", "澳门特别行政区"),
    "South China",
    province_name %in% c("四川省", "贵州省", "云南省", "西藏自治区", "重庆市"),
    "Southwest China",
    province_name %in% c("陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区"),
    "Northwest China"
  )
}

province <- st_read(province_shp_file) %>%
  st_make_valid() %>%
  mutate(region = assign_region(省))

province_plot <- province %>%
  rename(name = 省)

region_plot <- province_plot %>%
  group_by(region) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

region <- bind_rows(
  province %>%
    group_by(region) %>%
    summarise(geometry = st_union(geometry), .groups = "drop"),
  province %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    mutate(region = "China")
)

faod_raw <- qread(impact_file, nthreads = qn)

faod <- faod_raw %>%
  filter(year %in% 2001:2022) %>%
  mutate(
    Faod_results = map(Faod_results, function(adata) {
      
      temp <- adata %>%
        pivot_wider(
          names_from = faod_level,
          values_from = contains("%")
        ) %>%
        rast(type = "xyz", crs = "epsg:4326")
      
      fraction <- temp[["fraction"]]
      temp <- temp["%"]
      
      province_level <- exact_extract(
        temp,
        province,
        "weighted_mean",
        weights = fraction,
        progress = FALSE
      ) %>%
        bind_cols(province, .) %>%
        st_drop_geometry() %>%
        pivot_longer(
          cols = contains("%"),
          names_prefix = "weighted_mean.",
          names_to = c("name", "faod_level"),
          names_sep = "_"
        ) %>%
        pivot_wider()
      
      region_level <- exact_extract(
        temp,
        region,
        "weighted_mean",
        weights = fraction,
        progress = FALSE
      ) %>%
        bind_cols(region, .) %>%
        st_drop_geometry() %>%
        pivot_longer(
          cols = contains("%"),
          names_prefix = "weighted_mean.",
          names_to = c("name", "faod_level"),
          names_sep = "_"
        ) %>%
        pivot_wider()
      
      list(
        province_level = province_level,
        region_level = region_level
      )
    })
  )

saveRDS(faod, summary_file)

province_results <- faod %>%
  unnest_wider(Faod_results) %>%
  select(crop, year, crop_parent, province_level) %>%
  unnest(province_level)

region_results <- faod %>%
  unnest_wider(Faod_results) %>%
  select(crop, year, crop_parent, region_level) %>%
  unnest(region_level)

province_results_avg <- province_results %>%
  group_by(crop_parent, 省, faod_level) %>%
  summarise(
    across(c(`5%`, `50%`, `95%`), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

region_results_avg <- region_results %>%
  group_by(crop_parent, region, faod_level) %>%
  summarise(
    across(c(`5%`, `50%`, `95%`), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

write_xlsx(
  province_results_avg,
  file.path(table_dir, "faod_province_results_avg.xlsx")
)

write_xlsx(
  region_results_avg,
  file.path(table_dir, "faod_region_results_avg.xlsx")
)

plot_faod_map <- function(acrop, plot_legend, subtitle) {
  
  province_results_avg %>%
    filter(crop_parent == acrop, faod_level == 0.40) %>%
    inner_join(province_plot, by = c("省" = "name")) %>%
    ggplot() +
    geom_sf(
      aes(fill = `50%` * 100, geometry = geometry),
      color = NA,
      size = 0.2,
      show.legend = plot_legend
    ) +
    geom_sf(
      data = region_plot,
      fill = NA,
      color = "black",
      size = 0.8,
      show.legend = FALSE
    ) +
    scale_fill_gradientn(
      name = "Percentage change in SIF",
      limits = c(-10, 10),
      oob = scales::squish,
      colours = WrensBookshelf::WB_brewer(
        "BabyWrenAndTheGreatGift",
        direction = -1
      ),
      na.value = "grey90",
      guide = guide_colorbar(
        title.position = "top",
        barwidth = 15,
        label.position = "bottom"
      )
    ) +
    labs(tag = subtitle) +
    xlab(NULL) +
    ylab(NULL) +
    theme_ipsum_rc(
      grid = FALSE,
      axis = FALSE,
      base_size = 15
    ) +
    theme(
      legend.position = c(0.2, 0.1),
      legend.direction = "horizontal",
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.tag = element_text(size = 50),
      legend.title.align = 0.5
    )
}

plot_region_curve <- function(acrop, aregion) {
  
  region_results_avg %>%
    filter(crop_parent == acrop, region == aregion) %>%
    drop_na() %>%
    mutate(
      faod_level = as.numeric(faod_level),
      across(c(`5%`, `50%`, `95%`), ~ .x * 100)
    ) %>%
    filter(faod_level <= 0.5) %>%
    ggplot(aes(x = faod_level, y = `50%`)) +
    geom_vline(
      xintercept = 0.20,
      linetype = "dashed"
    ) +
    geom_ribbon(
      aes(ymin = `5%`, ymax = `95%`),
      alpha = 0.25,
      color = NA
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.6, color = "black") +
    labs(
      title = paste(acrop, aregion, sep = " - "),
      x = "Annual Faod level",
      y = "Percentage change in SIF"
    ) +
    theme_half_open(16, font_family = "Roboto Condensed") +
    background_grid() +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 18)
    )
}

crop_list <- unique(faod$crop_parent)

map_list <- pmap(
  list(
    crop_list,
    c(FALSE, TRUE, FALSE),
    letters[1:length(crop_list)]
  ),
  plot_faod_map
)

faod_map_patch <- wrap_plots(
  map_list,
  ncol = 2,
  byrow = FALSE
)

ggsave(
  file.path(figure_dir, "Figure_faod_impact_map.pdf"),
  faod_map_patch,
  width = 10,
  height = 12
)

ggsave(
  file.path(figure_dir, "Figure_faod_impact_map.tif"),
  faod_map_patch,
  width = 10,
  height = 12,
  dpi = 600,
  compression = "lzw"
)

region_curve_dir <- file.path(figure_dir, "region_curves")
dir.create(region_curve_dir, recursive = TRUE, showWarnings = FALSE)

region_curve_list <- expand_grid(
  crop_parent = crop_list,
  region = region_order
) %>%
  filter(region %in% unique(region_results_avg$region))

walk2(
  region_curve_list$crop_parent,
  region_curve_list$region,
  function(acrop, aregion) {
    
    p <- plot_region_curve(acrop, aregion)
    
    file_stub <- str_c(
      "Figure_faod_region_response_",
      str_replace_all(acrop, "[^A-Za-z0-9]+", "_"),
      "_",
      str_replace_all(aregion, "[^A-Za-z0-9]+", "_")
    )
    
    ggsave(
      file.path(region_curve_dir, str_c(file_stub, ".pdf")),
      p,
      width = 6,
      height = 4.5
    )
    
    ggsave(
      file.path(region_curve_dir, str_c(file_stub, ".tif")),
      p,
      width = 6,
      height = 4.5,
      dpi = 600,
      compression = "lzw"
    )
  }
)