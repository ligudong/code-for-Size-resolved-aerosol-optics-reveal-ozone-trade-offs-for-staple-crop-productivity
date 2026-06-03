source("R/00_load_packages.R")
source("R/00_helper_functions.R")

fpar_dir <- "data/raw/fpar"
mask_dir <- "data/derived/masks"
fpar_output_dir <- "data/derived/fpar"

dir.create(fpar_output_dir, recursive = TRUE, showWarnings = FALSE)

mask_path <- file.path(mask_dir, "mask_extraction.tif")
mask <- rast(mask_path)[[1]] %>% wrap()

plan(multisession, workers = 3)

fpar_data <- tibble(
  files_ = list.files(fpar_dir, pattern = "\\.tif$"),
  files  = list.files(fpar_dir, full.names = TRUE, pattern = "\\.tif$")
) %>%
  mutate(
    year  = substr(files_, 7, 10) %>% as.integer(),
    month = substr(files_, 11, 12) %>% as.integer(),
    date  = ymd(paste0(year, "-", month, "-01"))
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

fpar_data <- fpar_data %>%
  unnest(data) %>%
  trim_xy()

check_join(fpar_data)

fpar_data <- fpar_data %>%
  mutate(fpar = rowSums(select(., -(1:5)), na.rm = TRUE)) %>%
  select(1:5, fpar) %>%
  select(-files_)

saveRDS(
  fpar_data,
  file.path(fpar_output_dir, "fpar.rds")
)