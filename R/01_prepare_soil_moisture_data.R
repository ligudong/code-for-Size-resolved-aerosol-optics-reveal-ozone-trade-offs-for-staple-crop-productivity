# Load packages and helper functions
source("R/00_load_packages.R")
source("R/00_helper_functions.R")

# Define directories
soil_raw_dir   <- "data/raw/soil_moisture"
mask_dir       <- "data/derived/masks"
soil_output_dir <- "data/derived/soil_moisture"
dir.create(soil_output_dir, recursive = TRUE, showWarnings = FALSE)

mask_path <- file.path(mask_dir, "mask_extraction.tif")
mask <- rast(mask_path)[[1]] %>% wrap()

plan(multisession, workers = 3)

# Load root soil moisture bins
root_files <- list.files(file.path(soil_raw_dir, "root"), full.names = TRUE)
sm_bins_root <- tibble(
  files_ = root_files,
  files = root_files
) %>%
  mutate(
    data = future_map(
      files,
      function(afile) {
        r <- rast(afile)
        resample(r, mask, method = "near") %>%
          mask(mask) %>%
          as.data.frame(xy = TRUE)
      },
      .progress = TRUE,
      .options = furrr_options(packages = c("terra", "tidyverse"))
    ),
    .keep = "unused"
  )

plan(sequential)

# Unnest and trim
sm_bins_root <- sm_bins_root %>%
  unnest(data) %>%
  trim_xy()

# Aggregate by month/year and generate root_1 ~ root_9
root <- pro_map(sm_bins_root$data, function(adf) {
  adf %>%
    pivot_wider(
      names_from = c(year, month),
      values_from = c(
        `[0–100]`, `(100–200]`, `(200–300]`,
        `(300–400]`, `(400–500]`, `(500–600]`,
        `(600–700]`, `(700–800]`, `(800–10000]`
      )
    ) %>%
    rast(type = "xyz", crs = crs(mask)) %>%
    resample(mask, method = "near") %>%
    mask(mask) %>%
    as.data.frame(xy = TRUE) %>%
    pivot_longer(
      -c(x, y),
      names_to = c("bins", "year", "month"),
      names_sep = "_",
      names_transform = list(year = as.integer, month = as.integer)
    ) %>%
    pivot_wider(names_from = bins) %>%
    trim_xy()
})

root <- bind_rows(root) %>%
  `names<-`(c("x", "y", "year", "month", str_c("root_", 1:9)))

check_join(root)
saveRDS(root, file.path(soil_output_dir, "root.rds"))

# Load surface soil moisture bins
surface_files <- list.files(file.path(soil_raw_dir, "surface"), full.names = TRUE)
sm_bins_surface <- tibble(
  files_ = surface_files,
  files = surface_files
) %>%
  mutate(
    data = future_map(
      files,
      function(afile) {
        r <- rast(afile)
        resample(r, mask, method = "near") %>%
          mask(mask) %>%
          as.data.frame(xy = TRUE)
      },
      .progress = TRUE,
      .options = furrr_options(packages = c("terra", "tidyverse"))
    ),
    .keep = "unused"
  )

plan(sequential)
sm_bins_surface <- sm_bins_surface %>%
  unnest(data) %>%
  trim_xy()

# Define humidity bins
humidity_columns <- c("[0–5]", "(5–10]", "(10–15]", "(15–20]", "(20–25]",
                      "(25–30]", "(30–35]", "(35–40]", "(40–1000]")

# Aggregate by month/year and generate surface_1 ~ surface_9
surface <- sm_bins_surface %>%
  group_by(year, month) %>%
  summarise(
    data = list(reduce(data, function(df1, df2) {
      joined_df <- left_join(df1, df2, by = c("x", "y", "year", "month"), suffix = c("", ".y"))
      joined_df <- joined_df %>%
        mutate(across(all_of(humidity_columns),
                      ~ coalesce(.x, 0) + coalesce(tryCatch(get(paste0(cur_column(), ".y")), error = function(e) NULL), 0),
                      .names = "{col}")) %>%
        rename_with(~ sub("\\.x$", "", .), ends_with(".x")) %>%
        select(-ends_with(".y"))
      return(joined_df)
    }))
  ) %>%
  ungroup()

# Process each month into raster long format
surface_long <- pro_map(surface$data, function(adf) {
  adf %>%
    pivot_wider(
      names_from = c(year, month),
      values_from = all_of(humidity_columns)
    ) %>%
    rast(type = "xyz", crs = crs(mask)) %>%
    resample(mask, method = "near") %>%
    mask(mask) %>%
    as.data.frame(xy = TRUE) %>%
    pivot_longer(
      -c(x, y),
      names_to = c("bins", "year", "month"),
      names_sep = "_",
      names_transform = list(year = as.integer, month = as.integer)
    ) %>%
    pivot_wider(names_from = bins) %>%
    trim_xy()
})

surface_long <- bind_rows(surface_long) %>%
  `names<-`(c("x", "y", "year", "month", str_c("surface_", 1:9)))

check_join(surface_long)
saveRDS(surface_long, file.path(soil_output_dir, "surface.rds"))

# Compute monthly mean soil moisture (optional summary)
# Can use root.rds and surface.rds as input
mean_soil_moisture <- bind_rows(root, surface_long) %>%
  group_by(year, month) %>%
  summarise(across(starts_with(c("root_", "surface_")), mean, na.rm = TRUE)) %>%
  ungroup()

saveRDS(mean_soil_moisture, file.path(soil_output_dir, "mean.rds"))