source("R/00_load_packages.R")
source("R/00_helper_functions.R")

precipitation_dir <- "data/raw/precipitation"
mask_dir <- "data/derived/masks"
precipitation_output_dir <- "data/derived/precipitation"

dir.create(
  precipitation_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

mask_path <- file.path(mask_dir, "mask_extraction.tif")
mask <- rast(mask_path)[[1]] %>% wrap()

plan(multisession, workers = 3)

precipitation <- tibble(
  files_ = list.files(
    precipitation_dir,
    pattern = "\\.tif$"
  ),
  files = list.files(
    precipitation_dir,
    full.names = TRUE,
    pattern = "\\.tif$"
  )
) %>%
  mutate(
    year = str_extract(files_, "\\d{4}") %>% as.integer(),
    month = str_extract(files_, "(?<=-)\\d{2}") %>% as.integer()
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
      .options = furrr_options(
        packages = c("terra", "tidyverse")
      )
    ),
    .keep = "unused"
  )

plan(sequential)

precipitation <- precipitation %>%
  unnest(data) %>%
  trim_xy() %>%
  rename(prep = total_precipitation_sum)

check_join(precipitation)

saveRDS(
  precipitation,
  file.path(
    precipitation_output_dir,
    "tidied.rds"
  )
)