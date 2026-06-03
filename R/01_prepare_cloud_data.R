source("R/00_load_packages.R")
source("R/00_helper_functions.R")

library(terra)
library(furrr)
library(tidyverse)
library(lubridate)

cloud_dir <- "data/raw/cloud_fraction"
mask_dir  <- "data/derived/masks"
cloud_output_dir <- "data/derived/cloud"

dir.create(cloud_output_dir, recursive = TRUE, showWarnings = FALSE)

mask_path <- file.path(mask_dir, "mask_extraction.tif")
mask <- rast(mask_path)[[1]] %>% wrap()

plan(multisession, workers = 3)

# Load cloud_fraction monthly data
cloud_fraction <- tibble(
  files_ = list.files(file.path(cloud_dir, "cloud_fraction"), pattern = "\\.tif$"),
  files  = list.files(file.path(cloud_dir, "cloud_fraction"), full.names = TRUE, pattern = "\\.tif$")
) %>%
  mutate(
    year  = str_extract(files_, "\\d{4}"),
    month = str_extract(files_, "(?<=-)\\d{2}"),
    date  = ymd(paste0(year, "-", month, "-01")),
    year  = year(date),
    month = month(date)
  ) %>%
  mutate(
    data = future_map(
      files,
      function(afile) {
        mask_local <- rast(mask_path)
        rast(afile) %>%
          resample(mask_local, method = "average") %>%
          mask(mask_local) %>%
          as.data.frame(xy = TRUE)
      },
      .progress = TRUE,
      .options = furrr_options(packages = c("terra", "tidyverse"))
    ),
    .keep = "unused"
  )

# Restore sequential plan
plan(sequential)

cloud_fraction <- cloud_fraction %>%
  unnest(data) %>%
  trim_xy()

# Load Cloud_optical_depth monthly data
plan(multisession, workers = 3)

cloud_optical_depth <- tibble(
  files_ = list.files(file.path(cloud_dir, "Cloud_optical_depth"), pattern = "\\.tif$"),
  files  = list.files(file.path(cloud_dir, "Cloud_optical_depth"), full.names = TRUE, pattern = "\\.tif$")
) %>%
  mutate(
    year  = str_extract(files_, "\\d{4}"),
    month = str_extract(files_, "(?<=-)\\d{2}"),
    date  = ymd(paste0(year, "-", month, "-01")),
    year  = year(date),
    month = month(date)
  ) %>%
  mutate(
    data = future_map(
      files,
      function(afile) {
        mask_local <- rast(mask_path)
        rast(afile) %>%
          resample(mask_local, method = "average") %>%
          mask(mask_local) %>%
          as.data.frame(xy = TRUE)
      },
      .progress = TRUE,
      .options = furrr_options(packages = c("terra", "tidyverse"))
    ),
    .keep = "unused"
  )

plan(sequential)

cloud_optical_depth <- cloud_optical_depth %>%
  unnest(data) %>%
  trim_xy()

check_join(cloud_optical_depth)

# Merge cloud_fraction and cloud_optical_depth
MODIS <- inner_join(
  cloud_fraction,
  cloud_optical_depth,
  by = c("x", "y", "year", "month")
) %>%
  mutate(
    cloud_fraction = (cloud_fraction - mean(cloud_fraction, na.rm = TRUE)) / sd(cloud_fraction, na.rm = TRUE),
    cloud = cloud_fraction * cld_opd_acha
  ) %>%
  select(-cloud_fraction, -cld_opd_acha, -files_.x, -files_.y)

saveRDS(MODIS, file.path(cloud_output_dir, "tidied.rds"))