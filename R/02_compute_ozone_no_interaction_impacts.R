# File: R/02_compute_ozone_no_interaction_impacts.R

source("R/00_load_packages.R")
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

library(terra)
library(dplyr)
library(purrr)
library(qs)
library(matrixStats)
library(furrr)
library(tidyr)

data_dir <- "data/derived/crop/maize"
output_dir <- "data/derived/impacts"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load panel data
fnl_data <- get_data() %>%
  drop_na() %>%
  mutate(crop_parent = ifelse(crop %in% c("Rice(LR)","Rice(SR&ER)"), "Rice", crop))

# Load bootstrap coefficients
f1 <- read_rds(file.path(data_dir, "boots_f3.rds")) %>%
  mutate(crop = crop_parent) %>%
  bind_rows((.) %>% filter(crop=="Rice") %>% mutate(crop="Rice(LR)")) %>%
  mutate(crop = fifelse(crop=="Rice","Rice(SR&ER)",crop))

# Load non-interaction ozone raster data
ozone <- qread(file.path(output_dir,"impacts_ozone_AOT40_nojiaohu.qs"), nthreads = qn)

# Filter years and compute province/region weighted averages
ozone <- ozone %>%
  filter(year %in% 2001:2022) %>%
  mutate(
    ozone_results = map(ozone_results, function(adata){
      temp <- adata %>%
        pivot_wider(names_from = peak_level, values_from = contains("%")) %>%
        rast(type="xyz", crs="epsg:4326")
      
      fraction <- temp[["fraction"]]
      temp <- temp["%"]
      
      # Province-level weighted mean
      province_level <- exact_extract(temp, province, "weighted_mean", weights=fraction, progress=FALSE) %>%
        bind_cols(province, .) %>%
        st_drop_geometry() %>%
        pivot_longer(contains("%"), names_prefix="weighted_mean.", names_to=c("name","peak_level"), names_sep="_") %>%
        pivot_wider()
      
      # Region-level weighted mean
      region_level <- exact_extract(temp, region, "weighted_mean", weights=fraction, progress=FALSE) %>%
        bind_cols(region, .) %>%
        st_drop_geometry() %>%
        pivot_longer(contains("%"), names_prefix="weighted_mean.", names_to=c("name","peak_level"), names_sep="_") %>%
        pivot_wider()
      
      lst(province_level, region_level)
    }, .progress=TRUE)
  )

# Save results
saveRDS(ozone, file.path(output_dir,"impacts_ozone_summarised_nojiaohu.rds"))