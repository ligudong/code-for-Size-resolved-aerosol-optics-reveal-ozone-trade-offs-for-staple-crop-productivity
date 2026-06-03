source("R/00_load_packages.R")
source("R/00_helper_functions.R")

phenology_dir <- "data/raw/crop_phenology"
calendar_output_dir <- "data/derived/calendar"

dir.create(calendar_output_dir, showWarnings = FALSE, recursive = TRUE)

reference_raster <- rast(
  file.path(phenology_dir, "CHN_Maize_HE_2019.tif")
)

reference_extent <- ext(reference_raster)

phenology_stages <- tibble(
  file_pattern = c(
    "Maize_HE", "Maize_MA", "Maize_V3",
    "Rice\\(LR\\)_HE", "Rice\\(LR\\)_MA", "Rice\\(LR\\)_TR",
    "Rice\\(SR&ER\\)_HE", "Rice\\(SR&ER\\)_MA", "Rice\\(SR&ER\\)_TR",
    "Wheat_GR&EM", "Wheat_HE", "Wheat_MA"
  ),
  stage_name = c(
    "Maize_HE", "Maize_MA", "Maize_V3",
    "Rice(LR)_HE", "Rice(LR)_MA", "Rice(LR)_TR",
    "Rice(SR&ER)_HE", "Rice(SR&ER)_MA", "Rice(SR&ER)_TR",
    "Wheat_GR&EM", "Wheat_HE", "Wheat_MA"
  )
) %>%
  mutate(
    stage_name = str_replace_all(stage_name, "V3|TR", "GR&EM"),
    data = map2(file_pattern, stage_name, function(pattern, stage) {
      
      phenology_raster <- list.files(
        phenology_dir,
        pattern = pattern,
        full.names = TRUE
      ) %>%
        map(rast) %>%
        map(~ extend(.x, reference_extent)) %>%
        rast()
      
      names(phenology_raster) <- 2000:2019
      
      phenology_raster <- project(
        phenology_raster,
        "epsg:4326",
        method = "near"
      )
      
      target_grid <- rast(resolution = 0.5) %>%
        crop(phenology_raster, snap = "out")
      
      phenology_raster <- resample(
        phenology_raster,
        target_grid,
        method = "average"
      )
      
      phenology_data <- phenology_raster %>%
        as.data.frame(xy = TRUE, na.rm = FALSE) %>%
        pivot_longer(
          cols = -c(x, y),
          names_to = "year",
          values_to = "value",
          values_drop_na = TRUE,
          names_transform = list(year = as.integer)
        ) %>%
        mutate(
          date = round(value) %>%
            as.character() %>%
            str_c(year, "-", .) %>%
            parse_date_time("Y-j"),
          month = month(date),
          day = day(date)
        )
      
      if (str_detect(stage, "GR&EM")) {
        phenology_data <- phenology_data %>%
          mutate(month_fnl = fifelse(day >= 15, month + 1, month))
      }
      
      if (str_detect(stage, "MA")) {
        phenology_data <- phenology_data %>%
          mutate(month_fnl = fifelse(day >= 15, month, month - 1))
      }
      
      if (str_detect(stage, "HE")) {
        phenology_data <- phenology_data %>%
          mutate(month_fnl = month)
      }
      
      phenology_data %>%
        select(x, y, year, month_fnl)
    })
  ) %>%
  select(-file_pattern)

crop_calendar <- phenology_stages %>%
  separate(stage_name, into = c("crop", "stage"), sep = "_") %>%
  unnest(data) %>%
  pivot_wider(
    names_from = "stage",
    values_from = "month_fnl"
  ) %>%
  mutate(month = map2(`GR&EM`, MA, seq)) %>%
  unnest(month) %>%
  filter(MA >= `GR&EM`)

saveRDS(
  crop_calendar,
  file.path(calendar_output_dir, "crop_calendar_0.5deg.rds")
)