source("R/00_load_packages.R")
source("R/00_helper_functions.R")

pm25_dir <- "data/raw/pm25"
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

PM25 <- tibble(
  files_ = list.files(
    pm25_dir,
    pattern = "\\.tif$"
  ),
  files = list.files(
    pm25_dir,
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

PM25 <- PM25 %>%
  unnest(data) %>%
  trim_xy()

check_join(PM25)

PM25 <- PM25 %>%
  mutate(
    PM25 = rowSums(
      select(., -(1:5)),
      na.rm = TRUE
    )
  ) %>%
  select(
    x,
    y,
    year,
    month,
    PM25
  )

saveRDS(
  PM25,
  file.path(
    pm_output_dir,
    "PM25.rds"
  )
)