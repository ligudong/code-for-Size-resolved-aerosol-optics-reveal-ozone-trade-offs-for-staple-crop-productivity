# Load packages and helper functions
source("R/00_load_packages.R")
source("R/00_helper_functions.R")

library(terra)
library(furrr)
library(tidyverse)
library(lubridate)

# Directories
mask_dir <- "data/derived/masks"
faod_dir <- "data/raw/faod_monthly_mean"
faod_output_dir <- "data/derived/aerosol"
dir.create(faod_output_dir, recursive = TRUE, showWarnings = FALSE)

# Load maize mask
mask <- rast(file.path(mask_dir, "mask_extraction.tif"))[[1]] %>% wrap()

# Set parallel processing
plan(multisession, workers = 3)

# Load Faod monthly files and extract year/month
Faod <- tibble(
  files_ = list.files(faod_dir, pattern = "\\.tif$"),
  files  = list.files(faod_dir, full.names = TRUE, pattern = "\\.tif$")
) %>%
  mutate(
    year  = substr(files_, 1, 4),
    month = substr(files_, 5, 6),
    date  = ymd(paste0(year, "-", month, "-01")),
    year  = year(date),
    month = month(date)
  ) %>%
  mutate(
    data = future_map(
      files,
      function(afile) {
        r <- rast(afile)
        # Resample to 0.5° and apply mask
        resample(r, rast(mask), method = "average") %>%
          mask(rast(mask)) %>%
          as.data.frame(xy = TRUE)
      },
      .progress = TRUE,
      .options = furrr_options(packages = c("terra", "tidyverse"))
    ),
    .keep = "unused"
  )

# Restore sequential plan
plan(sequential)

# Unnest data and round coordinates
Faod <- Faod %>%
  unnest(data) %>%
  trim_xy()

# Check grid overlap
check_join(Faod)

# Sum all Faod columns (after first five: x, y, year, month, date)
Faod <- Faod %>%
  mutate(Faod = rowSums(select(., -(1:5)), na.rm = TRUE)) %>%
  select(1:5, Faod)

# Remove temporary columns
Faod <- Faod %>% select(-files_)

# Save processed Faod
saveRDS(Faod, file.path(faod_output_dir, "Faod.rds"))