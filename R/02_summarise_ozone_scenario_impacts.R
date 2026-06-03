source("R/00_load_packages.R")
source("R/00_helper_functions.R")

province_shp_file <- "data/raw/shp/province.shp"
impact_file_base <- "data/impacts_ozone_AOT40.qs"
impact_file_nointeraction <- "data/impacts_ozone_AOT40_nojiaohu.qs"

output_dir <- "data/derived/ozone_summary_tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Assign region helper
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

province <- st_read(province_shp_file) %>%
  st_make_valid() %>%
  mutate(region = assign_region(省))

region <- bind_rows(
  province %>% group_by(region) %>% summarise(geometry = st_union(geometry), .groups="drop"),
  province %>% summarise(geometry = st_union(geometry), .groups="drop") %>% mutate(region="China")
)

process_ozone_file <- function(file_path, output_rds) {
  ozone <- qread(file_path, nthreads=qn) %>%
    filter(year %in% 2001:2022) %>%
    mutate(
      ozone_results = map(ozone_results, function(adata) {
        temp <- adata %>%
          pivot_wider(names_from = peak_level, values_from = contains("%")) %>%
          rast(type="xyz", crs="epsg:4326")
        fraction <- temp[["fraction"]]
        temp <- temp["%"]
        
        # Province-level summary
        province_level <- exact_extract(temp, province, "weighted_mean", weights=fraction, progress=FALSE) %>%
          bind_cols(province, .) %>%
          st_drop_geometry() %>%
          pivot_longer(contains("%"), names_prefix="weighted_mean.", names_to=c("name","peak_level"), names_sep="_") %>%
          pivot_wider()
        
        # Region-level summary
        region_level <- exact_extract(temp, region, "weighted_mean", weights=fraction, progress=FALSE) %>%
          bind_cols(region, .) %>%
          st_drop_geometry() %>%
          pivot_longer(contains("%"), names_prefix="weighted_mean.", names_to=c("name","peak_level"), names_sep="_") %>%
          pivot_wider()
        
        lst(province_level, region_level)
      }, .progress=TRUE)
    )
  
  saveRDS(ozone, output_rds)
  return(ozone)
}

# Process base interaction O3
ozone_base <- process_ozone_file(impact_file_base, file.path(output_dir,"impacts_ozone_summarised_base.rds"))

# Process no interaction O3
ozone_nointer <- process_ozone_file(impact_file_nointeraction, file.path(output_dir,"impacts_ozone_summarised_no_interaction.rds"))