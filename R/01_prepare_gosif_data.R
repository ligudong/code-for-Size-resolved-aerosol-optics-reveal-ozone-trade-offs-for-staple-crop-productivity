source("R/00_load_packages.R")
source("R/00_helper_functions.R")

gosif_dir <- "data/raw/gosif"
mask_dir <- "data/derived/masks"
sif_output_dir <- "data/derived/sif"

dir.create(
  sif_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

mask <- rast(
  file.path(mask_dir, "mask_extraction.tif")
)[[1]] %>%
  wrap()

plan(multisession, workers = 3)

GOSIF <- tibble(
  files_ = list.files(
    gosif_dir,
    pattern = "\\.tif$"
  ),
  files = list.files(
    gosif_dir,
    pattern = "\\.tif$",
    full.names = TRUE
  )
) %>%
  mutate(
    year = str_extract(
      files_,
      "(?<=GOSIF_)\\d{4}"
    ),
    month = str_extract(
      files_,
      "(?<=\\.M)\\d{2}"
    )
  ) %>%
  mutate(
    date = ymd(
      paste0(
        year,
        "-",
        month,
        "-01"
      )
    )
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
        
        r <- rast(afile)
        
        resample(
          r,
          rast(mask),
          method = "average"
        ) %>%
          mask(rast(mask)) %>%
          as.data.frame(xy = TRUE)
      },
      .progress = TRUE,
      .options = furrr_options(
        packages = c(
          "terra",
          "tidyverse"
        )
      )
    ),
    .keep = "unused"
  )

plan(sequential)

GOSIF <- GOSIF %>%
  unnest(data) %>%
  trim_xy()

check_join(GOSIF)

GOSIF <- GOSIF %>%
  mutate(
    GOSIF = rowSums(
      select(
        .,
        starts_with("GOSIF_")
      ),
      na.rm = TRUE
    )
  ) %>%
  select(
    -starts_with("GOSIF_")
  ) %>%
  select(
    -files_
  )

saveRDS(
  GOSIF,
  file.path(
    sif_output_dir,
    "GOSIF.rds"
  )
)