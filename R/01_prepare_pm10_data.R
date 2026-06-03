source("R/00_load_packages.R")
source("R/00_helper_functions.R")

pm10_dir <- "data/raw/pm10"
mask_dir <- "data/derived/masks"
pm_output_dir <- "data/derived/pm"

dir.create(
  pm_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

mask_path <- file.path(mask_dir, "mask_extraction.tif")

mask <- rast(mask_path) %>%
  sum(na.rm = TRUE) %>%
  wrap()

plan(multisession, workers = 3)

PM10 <- tibble(
  files_ = list.files(
    pm10_dir,
    pattern = "\\.tif$"
  ),
  files = list.files(
    pm10_dir,
    full.names = TRUE,
    pattern = "\\.tif$"
  )
) %>%
  mutate(
    year = str_extract(files_, "\\d{4}") %>% as.integer(),
    month = str_extract(files_, "(?<!\\d)(0[1-9]|1[0-2])(?!\\d)") %>% as.integer()
  ) %>%
  mutate(
    data = future_map(
      files,
      function(afile) {
        mask_local <- rast(mask_path) %>%
          sum(na.rm = TRUE)
        
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

PM10 <- PM10 %>%
  unnest(data) %>%
  trim_xy()

check_join(PM10)

PM10 <- PM10 %>%
  mutate(
    PM10 = rowSums(
      select(., -(1:5)),
      na.rm = TRUE
    )
  ) %>%
  select(
    x,
    y,
    year,
    month,
    PM10
  )

saveRDS(
  PM10,
  file.path(
    pm_output_dir,
    "PM10.rds"
  )
)