source("R/00_load_packages.R")
source("R/00_helper_functions.R")

rtsif_dir <- "data/raw/rtsif"
mask_dir <- "data/derived/masks"
sif_output_dir <- "data/derived/sif"

dir.create(sif_output_dir, recursive = TRUE, showWarnings = FALSE)

mask_path <- file.path(mask_dir, "mask_extraction.tif")
mask <- rast(mask_path) %>% wrap()

plan(multisession, workers = 3)

rtSIF <- tibble(
  files_ = list.files(rtsif_dir, pattern = "\\.tif$"),
  files  = list.files(rtsif_dir, full.names = TRUE, pattern = "\\.tif$")
) %>%
  mutate(
    year = str_extract(files_, "\\d{4}"),
    month = str_extract(files_, "(?<=-)[0-9]{2}(?=-mean\\.tif)"),
    date = ymd(paste0(year, "-", month, "-01"))
  ) %>%
  mutate(year = year(date), .keep = "unused") %>%
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

rtSIF <- rtSIF %>%
  unnest() %>%
  trim_xy()

# 合并所有图层为一列
rtSIF <- rtSIF %>%
  mutate(rtsif = coalesce(lyr1, lyr2, lyr3)) %>%
  select(x, y, year, month, rtsif) %>%
  distinct(x, y, year, month, rtsif, .keep_all = TRUE)

saveRDS(rtSIF, file.path(sif_output_dir, "RTSIF.rds"))