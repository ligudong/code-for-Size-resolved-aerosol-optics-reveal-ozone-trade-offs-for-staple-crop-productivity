# Load packages and helper functions
source("R/00_load_packages.R")
source("R/00_helper_functions.R")

# Directories
mask_dir <- "data/derived/masks"
caod_dir <- "data/raw/caod_monthly_mean"
caod_output_dir <- "data/derived/aerosol"
dir.create(caod_output_dir, recursive = TRUE, showWarnings = FALSE)

# Load maize mask
mask <- rast(file.path(mask_dir, "mask_extraction.tif"))[[1]] %>% wrap()

# Set parallel processing
plan(multisession, workers = 3)

# Load monthly Caod files
Caod <- tibble(
  files_ = list.files(caod_dir, pattern = "\\.tif$"),
  files  = list.files(caod_dir, full.names = TRUE, pattern = "\\.tif$")
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

# Unnest data and trim coordinates
Caod <- Caod %>%
  unnest(data) %>%
  trim_xy()

# Sum all Caod columns after the first 5 (x, y, year, month, date)
Caod <- Caod %>%
  mutate(Caod = rowSums(select(., -(1:5)), na.rm = TRUE)) %>%
  select(2:5, Caod)  # Keep x, y, year, month, Caod

# Save processed Caod
saveRDS(Caod, file.path(caod_output_dir, "Caod.rds"))