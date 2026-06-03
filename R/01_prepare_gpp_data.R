source("R/00_load_packages.R")
source("R/00_helper_functions.R")

gpp_dir <- "data/raw/gpp"
mask_dir <- "data/derived/masks"
gpp_output_dir <- "data/derived/gpp"

dir.create(
  gpp_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

mask_path <- file.path(mask_dir, "mask_extraction.tif")
mask <- rast(mask_path) %>% wrap()

plan(multisession, workers = 3)

GPP <- tibble(
  files_ = list.files(
    gpp_dir,
    pattern = "\\.tif$"
  ),
  files = list.files(
    gpp_dir,
    full.names = TRUE,
    pattern = "\\.tif$"
  )
) %>%
  mutate(
    year = str_extract(files_, "\\d{4}"),
    date = ymd(paste0(year, "-01-01")),
    year = year(date),
    month = 1
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

GPP <- GPP %>%
  unnest(data) %>%
  trim_xy()

GPP <- GPP %>%
  mutate(
    GPP = coalesce(lyr1, lyr2, lyr3)
  ) %>%
  select(
    x,
    y,
    year,
    GPP
  ) %>%
  distinct(
    x,
    y,
    year,
    GPP,
    .keep_all = TRUE
  )

check_join(GPP)

saveRDS(
  GPP,
  file.path(
    gpp_output_dir,
    "GPP.rds"
  )
)