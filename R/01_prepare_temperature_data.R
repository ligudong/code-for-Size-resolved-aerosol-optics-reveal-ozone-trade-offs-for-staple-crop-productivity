source("R/00_load_packages.R")
source("R/00_helper_functions.R")

library(terra)
library(furrr)
library(tidyverse)
library(lubridate)
library(qs)

# Directories
mask_dir <- "data/derived/masks"
temp_dir <- "data/raw/temp"
temp_output_dir <- "data/derived/temp"

dir.create(temp_output_dir, recursive = TRUE, showWarnings = FALSE)

mask_path <- file.path(mask_dir, "mask_extraction.tif")
mask <- rast(mask_path) %>% wrap()

plan(multisession, workers = 3)

# Load temperature monthly files
temp_files <- list.files(temp_dir, pattern = "\\.tif$", full.names = TRUE)
temp_names <- list.files(temp_dir, pattern = "\\.tif$")

temp <- tibble(
  files_ = temp_names,
  files  = temp_files
) %>%
  mutate(
    year = str_extract(files_, "\\d{4}"),
    date = str_extract(files_, "\\d{4}-\\d{2}-\\d{2}") %>% ymd()
  ) %>%
  mutate(
    year = year(date),
    month = month(date),
    .keep = "unused"
  ) %>%
  mutate(
    data = future_map(
      files,
      function(afile) {
        mask_local <- rast(mask_path)
        r <- rast(afile)
        resample(r, mask_local, method = "average") %>%
          mask(mask_local) %>%
          as.data.frame(xy = TRUE)
      },
      .progress = TRUE,
      .options = furrr_options(packages = c("terra", "tidyverse"))
    ),
    .keep = "unused"
  )

plan(sequential)

# Unnest and trim
temp <- temp %>%
  unnest(data) %>%
  trim_xy() %>%
  rename(temp = last_col())

# Save raw temp data
saveRDS(temp, file.path(temp_output_dir, "tmax.rds"))

# Convert Kelvin to Celsius
temp <- temp %>%
  mutate(temp_celsius = temp - 273.15)

# Define bins
bins <- c(-Inf, seq(0, 40, by = 1), Inf)
bin_labels <- c("<0", as.character(0:39), ">40")

# Bin counts per month
temperature_summary <- temp %>%
  mutate(
    bin = cut(temp_celsius, breaks = bins, labels = bin_labels, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(year, month, x, y, bin) %>%
  summarise(freq = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = bin,
    values_from = freq,
    values_fill = 0
  )

# Rename columns to temp_1 ~ temp_42
temperature_summary <- temperature_summary %>%
  rename_with(~ c("temp_1", paste0("temp_", 2:41), "temp_42"), .cols = everything())

# Save processed binned temperature
qsave(temperature_summary, file.path(temp_output_dir, "temp_bins.qs"), nthreads = 5)