source("R/00_load_packages.R")
source("R/00_helper_functions.R")

fertilizer_dir <- "data/raw/fertilizer"
mask_dir <- "data/derived/masks"
fertilizer_output_dir <- "data/derived/fertilizer"

dir.create(
  fertilizer_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

mask_path <- file.path(mask_dir, "mask_extraction.tif")
mask <- rast(mask_path)[[1]] %>% wrap()

plan(multisession, workers = 3)

fertilizer_data <- tibble(
  files_ = list.files(
    fertilizer_dir,
    pattern = "\\.tif$"
  ),
  files = list.files(
    fertilizer_dir,
    full.names = TRUE,
    pattern = "\\.tif$"
  )
) %>%
  mutate(
    type = str_extract(files_, "^[A-Za-z]+"),
    position = str_extract(files_, "(deep|surface)"),
    year = str_extract(files_, "\\d{4}") %>% as.integer()
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

fertilizer_data <- fertilizer_data %>%
  unnest(data) %>%
  trim_xy()

fertilizer_data <- fertilizer_data %>%
  mutate(
    Total = rowSums(
      select(., 7:ncol(.)),
      na.rm = TRUE
    )
  ) %>%
  select(
    x,
    y,
    type,
    position,
    year,
    Total
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = Total,
    values_fill = 0
  )

fertilizer_deep <- fertilizer_data %>%
  filter(position == "deep")

fertilizer_surface <- fertilizer_data %>%
  filter(position == "surface")

saveRDS(
  fertilizer_deep,
  file.path(fertilizer_output_dir, "fertilizer_deep.rds")
)

saveRDS(
  fertilizer_surface,
  file.path(fertilizer_output_dir, "fertilizer_surface.rds")
)