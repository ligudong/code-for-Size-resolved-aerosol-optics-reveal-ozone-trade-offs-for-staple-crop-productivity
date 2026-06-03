source("R/00_load_packages.R")
source("R/00_helper_functions.R")

library(terra)
library(tidyverse)
library(dtplyr)
library(sf)
library(qs)
library(lubridate)

# Define paths
mask_dir        <- "data/derived/masks"
calendar_file   <- "data/derived/calendar/tidied.rds"
fraction_dir    <- "data/derived/masks"
sif_dir         <- "data/derived/sif"
temp_dir        <- "data/derived/temp"
soil_dir        <- "data/derived/soil_moisture"
cloud_dir       <- "data/derived/cloud"
aerosol_dir     <- "data/derived/aerosol"
ozone_dir       <- "data/derived/ozone"
pm_dir          <- "data/derived/pm"
vpd_dir         <- "data/derived/vpd"
fertilizer_dir  <- "data/derived/fertilizer"
irrigation_dir  <- "data/derived/irrigation"
yield_dir       <- "data/derived/crop_production"
co2_file        <- "data/derived/co2.rds"
shp_file_path   <- "data/raw/shp/county.shp"

# Load panel data modules
cldr     <- read_rds(calendar_file) %>% filter((MA - `GR&EM`) >= 2) %>% trim_xy()
fraction <- map(c("Maize", "Rice", "Rice", "Wheat"), ~{
  afile <- file.path(fraction_dir, str_c("mask_", .x, ".tif"))
  rast(afile) %>%
    as.data.frame(xy = TRUE, na.rm = FALSE) %>%
    pivot_longer(-c(x, y),
                 names_to = "year", names_transform = list(year = as.integer),
                 values_drop_na = TRUE,
                 values_to = "fraction")
}) %>% bind_rows() %>% trim_xy()

GOSIF  <- read_rds(file.path(sif_dir, "GOSIF.rds")) %>% mutate(GOSIF = GOSIF*0.0001)
CSIF   <- read_rds(file.path(sif_dir, "CSIF.rds"))
RTSIF  <- read_rds(file.path(sif_dir, "RTSIF.rds"))

tmax    <- read_rds(file.path(temp_dir, "tmax.rds")) %>% mutate(maxtmp = temp - 273.15) %>% select(-files_, -temp)
temp_bins <- read_qs(file.path(temp_dir, "temp_bins.qs")) %>% rename_with(~str_replace(.x, "^temp_", "bins"), starts_with("temp"))
prep    <- read_rds(file.path(temp_dir, "tidied.rds")) %>% select(-files_)
surface <- read_rds(file.path(soil_dir, "surface.rds"))
root    <- read_rds(file.path(soil_dir, "root.rds"))
sm_root <- read_rds(file.path(soil_dir, "root_data.rds")) %>% select(-files_) %>% group_by(x, y, year, month) %>% summarise(mean_root = mean(root, na.rm = TRUE))
sm_surface <- read_rds(file.path(soil_dir, "surface_data.rds")) %>% group_by(x, y, year, month) %>% summarise(mean_surface = mean(surface, na.rm = TRUE))
cloud   <- read_rds(file.path(cloud_dir, "tidied.rds"))
Caod    <- read_rds(file.path(aerosol_dir, "Caod.rds"))
Faod    <- read_rds(file.path(aerosol_dir, "Faod.rds"))
ozone   <- read_rds(file.path(ozone_dir, "tidied.rds"))

# Load additional environmental data: PM2.5, PM10, VPD
PM25    <- readRDS(file.path(pm_dir, "PM25.rds"))
PM10    <- readRDS(file.path(pm_dir, "PM10.rds"))
VPD     <- readRDS(file.path(vpd_dir, "tidied.rds"))

co2     <- read_rds(co2_file)
irrigation <- read_rds(file.path(irrigation_dir, "tidied.rds")) %>% fgroup_by(x, y) %>% fsummarise(irg_fraction = fmean(irg_fraction)) %>% fungroup()

fertilizer_deep    <- read_rds(file.path(fertilizer_dir, "fertilizer_deep.rds")) %>% select(-position) %>% rename(D1=AA,D2=AN,D3=AP,D4=AS,D5=CAN,D6=CR,D7=MA,D8=NK,D9=NPK,D10=NS,D11=ONP,D12=ONS,D13=Urea)
fertilizer_surface <- read_rds(file.path(fertilizer_dir, "fertilizer_surface.rds")) %>% select(-position) %>% rename(d1=AA,d2=AN,d3=AP,d4=AS,d5=CAN,d6=CR,d7=MA,d8=NK,d9=NPK,d10=NS,d11=ONP,d12=ONS,d13=Urea)
yield <- read_rds(file.path(yield_dir, "crop_production.rds"))

# Ensure year/month are integers
list(cldr, fraction, GOSIF, temp_bins, tmax, prep, surface, root, cloud, Faod, Caod, sm_root, sm_surface, co2, PM25, PM10, VPD) %>%
  walk(~ mutate(.x, year = as.integer(year)))
list(GOSIF, temp_bins, tmax, prep, surface, root, cloud, Faod, Caod, sm_root, sm_surface, PM25, PM10, VPD) %>%
  walk(~ mutate(.x, month = as.integer(month)))

# Merge all data
data <- cldr %>%
  inner_join(fraction) %>%
  left_join(GOSIF) %>%
  left_join(CSIF) %>%
  left_join(RTSIF) %>%
  left_join(temp_bins) %>%
  left_join(tmax) %>%
  left_join(prep) %>%
  left_join(surface) %>%
  left_join(root) %>%
  left_join(sm_root) %>%
  left_join(sm_surface) %>%
  left_join(cloud) %>%
  left_join(fertilizer_deep) %>%
  left_join(fertilizer_surface) %>%
  left_join(Faod) %>%
  left_join(Caod) %>%
  left_join(PM25) %>%
  left_join(PM10) %>%
  left_join(VPD)

# Load county shapefile and generate xy info
shp_file <- st_read(file.path("data/raw/shp", "county.shp")) %>% transmute(
  county = NAME, county_code = PAC,
  city = 市, city_code = 市代码,
  province = 省, province_code = 省代码
)

xy_info <- data %>%
  distinct(x, y) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(shp_file), remove = FALSE) %>%
  st_intersection(shp_file) %>%
  st_drop_geometry()

# Ozone summary
ozone_s <- map(1:7, function(anum) {
  cldr %>%
    filter((MA - month) < anum) %>%
    inner_join(fraction) %>%
    left_join(ozone) %>%
    lazy_dt() %>%
    group_by(x, y, year, crop) %>%
    summarise(
      O3 = sum(O3, na.rm = TRUE),
      across(c(W126, AOT40), sum, na.rm = TRUE)
    ) %>%
    as_tibble() %>%
    rename_with(~ str_c(.x, "_", anum), c(O3, W126, AOT40))
}) %>%
  reduce(full_join)

# Aggregate variables to generate final fnl_data
fnl_data <- data %>%
  lazy_dt() %>%
  group_by(x, y, year, crop) %>%
  summarise(
    across(c(GOSIF), max, .names = "{.col}_peak"),
    across(c(GOSIF), sum, .names = "{.col}_sum"),
    across(c(cloud, Faod, Caod, HE, MA, `GR&EM`, maxtmp, mean_root, mean_surface, GOSIF, fraction, PM25, PM10, VPD), mean),
    across(c(prep, starts_with("bin"), starts_with("step"), starts_with("D"), starts_with("d"), starts_with("surface_"), starts_with("root_")), sum)
  ) %>%
  mutate(x_y = str_c(x, "_", y)) %>%
  inner_join(xy_info) %>%
  inner_join(irrigation) %>%
  left_join(ozone_s) %>%
  as_tibble()

fnl_data <- fnl_data %>% filter(!is.na(GOSIF))
fnl_data[] <- lapply(fnl_data, function(x) ifelse(is.na(x), 0, x))
fnl_data <- fnl_data %>%
  rename_with(~ str_replace_all(., "bins", "bin"))

qsave(fnl_data, file.path("data/derived", "tidied.qs"), nthreads = 10)