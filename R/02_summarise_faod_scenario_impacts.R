source("R/00_load_packages.R")
source("R/00_helper_functions.R")

# Directories
province_shp_file <- "data/raw/shp/province.shp"
impact_file      <- "data/derived/impacts/impacts_Faod.qs"
output_dir       <- "data/derived/faod_summary_tables"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Assign regions
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

# Load province shapefile
province <- st_read(province_shp_file) %>%
  st_make_valid() %>%
  mutate(region = assign_region(省))

# Combine provinces into regions
region <- bind_rows(
  province %>%
    group_by(region) %>%
    summarise(geometry = st_union(geometry), .groups = "drop"),
  province %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    mutate(region = "China")
)

# Load Faod impact data
faod <- qread(impact_file, nthreads = qn) %>%
  filter(year %in% 2001:2022)

# Compute province-level and region-level weighted mean
faod_summary <- faod %>%
  mutate(
    Faod_results = map(Faod_results, function(adata) {
      temp <- adata %>%
        pivot_wider(names_from = faod_level, values_from = contains("%")) %>%
        rast(type = "xyz", crs = "epsg:4326")
      
      fraction <- temp[["fraction"]]
      temp <- temp["%"]
      
      region_level <- exact_extract(temp, region, "weighted_mean", weights = fraction, progress = FALSE) %>%
        bind_cols(region, .) %>%
        st_drop_geometry() %>%
        pivot_longer(
          cols = contains("%"),
          names_prefix = "weighted_mean.",
          names_to = c("name", "faod_level"),
          names_sep = "_"
        ) %>%
        pivot_wider()
      
      region_level
    })
  )

# Unnest and compute yearly mean for regions
region_year_table <- faod_summary %>%
  select(crop, crop_parent, year, Faod_results) %>%
  unnest(Faod_results) %>%
  mutate(faod_level = as.numeric(faod_level))

region_average_table <- region_year_table %>%
  group_by(crop_parent, region, faod_level) %>%
  summarise(
    across(c(`5%`, `50%`, `95%`), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

china_average_table <- region_average_table %>%
  filter(region == "China")

# Export Excel and RDS files
write_xlsx(
  list(
    region_year = region_year_table,
    region_average = region_average_table,
    china_average = china_average_table
  ),
  file.path(output_dir, "faod_summary_tables.xlsx")
)

saveRDS(region_year_table, file.path(output_dir, "faod_region_year_table.rds"))
saveRDS(region_average_table, file.path(output_dir, "faod_region_average_table.rds"))
saveRDS(china_average_table, file.path(output_dir, "faod_china_average_table.rds"))