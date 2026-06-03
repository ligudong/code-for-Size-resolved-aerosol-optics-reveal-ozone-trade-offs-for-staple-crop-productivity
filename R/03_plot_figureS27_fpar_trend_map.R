source("R/00_load_packages.R")
source("R/00_helper_functions.R")

data_path <- file.path("data/derived", "tidied.qs")
figure_dir <- "figures/supplementary/figureS27"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

region_order <- c(
  "North China", "Northeast China", "East China", "South China",
  "Central China", "Southwest China", "Northwest China"
)

province_plot <- st_read(
  "https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json",
  quiet = TRUE
) %>%
  st_make_valid() %>%
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

region_plot <- province_plot %>%
  group_by(region) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

city <- province_plot %>%
  st_drop_geometry() %>%
  filter(adcode != "100000_JD") %>%
  transmute(
    adcode = fifelse(
      adcode %in% c("110000", "120000", "310000", "500000", "460000", "710000", "810000", "820000"),
      adcode,
      str_c(adcode, "_full")
    ),
    data = map(
      adcode,
      ~ st_read(
        str_c("https://geo.datav.aliyun.com/areas_v3/bound/", .x, ".json"),
        quiet = TRUE
      )
    )
  ) %>%
  as_tibble()

city <- bind_rows(
  city %>% filter(str_detect(adcode, "full")) %>% unnest(data),
  city %>% filter(!str_detect(adcode, "full")) %>% unnest(data)
) %>%
  st_as_sf(sf_column_name = "geometry") %>%
  st_make_valid()

data <- qread(data_path, nthreads = qn) %>%
  mutate(
    crop = fifelse(str_detect(crop, "Rice"), "Rice", crop)
  )

yeartrend <- data %>%
  group_by(x, y, crop) %>%
  summarise(
    across(
      c(contains("GOSIF"), fraction, Faod, Caod, fpar),
      ~ tryCatch(
        coef(lm(. ~ year, na.action = na.exclude))[2],
        error = function(e) NA_real_
      ),
      .names = "trend_{.col}"
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = crop,
    values_from = starts_with("trend_")
  ) %>%
  as_tibble()

if (!"trend_fpar_Maize" %in% names(yeartrend)) {
  stop("trend_fpar_Maize column not found. Please check the input data.")
}

yeartrend_vect <- vect(
  yeartrend,
  geom = c("x", "y"),
  crs = "epsg:4326"
)

yeartrend_rast <- rasterize(
  yeartrend_vect,
  rast(yeartrend_vect, res = 0.25),
  field = "trend_fpar_Maize"
)

extracted <- exact_extract(
  yeartrend_rast,
  city,
  "mean"
)

extracted <- bind_cols(city, extracted) %>%
  rename(fpar_trend_Maize = mean) %>%
  mutate(
    area = units::drop_units(st_area(.)) / 1e6
  )

custom_colors <- c(
  "#20223e", "#3b1f46", "#7f4b89", "#b46db3", "#e3a5d6",
  "#f7f7f7", "white", "#f7f7f7",
  "#fcc893", "#feb424", "#fd8700"
)

color_values <- c(
  -0.02, -0.015, -0.01, -0.005, -0.002,
  -0.0001, 0, 0.0001,
  0.005, 0.01, 0.015
)

plot_s27 <- ggplot(extracted) +
  geom_sf(
    aes(fill = fpar_trend_Maize),
    color = NA
  ) +
  geom_sf(
    data = region_plot,
    fill = NA,
    color = "black",
    linewidth = 1.2
  ) +
  scale_fill_gradientn(
    name = "FPAR trend (Maize)",
    colours = custom_colors,
    values = scales::rescale(color_values, to = c(0, 1)),
    limits = c(-0.02, 0.015),
    oob = scales::squish,
    na.value = "grey90",
    guide = guide_colorbar(
      title.position = "top",
      barwidth = 10,
      label.position = "bottom"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.2, 0.2),
    legend.direction = "horizontal",
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

print(plot_s27)

ggsave(
  filename = file.path(figure_dir, "figureS27_fpar_trend_map_maize.tif"),
  plot = plot_s27,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw"
)