source("R/00_load_packages.R")
source("R/00_helper_functions.R")

production_dir <- "data/raw/crop_production"
production_output_dir <- "data/derived/crop_production"

dir.create(
  production_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

crops <- c("Maize", "Rice", "Wheat")

read_crop_production <- function(crop_name) {
  
  crop_path <- file.path(production_dir, crop_name)
  
  files <- list.files(
    path = crop_path,
    pattern = paste0("GGCP10_Production_.*_", crop_name, "\\.tif$"),
    full.names = TRUE
  )
  
  map_dfr(
    files,
    function(afile) {
      
      year_current <- str_extract(
        basename(afile),
        paste0("(?<=_)\\d{4}(?=_{0,1}", crop_name, ")")
      ) %>%
        as.integer()
      
      rast(afile) %>%
        as.data.frame(xy = TRUE) %>%
        rename(yield = last_col()) %>%
        mutate(
          year = year_current,
          crop = crop_name
        ) %>%
        select(
          crop,
          year,
          x,
          y,
          yield
        ) %>%
        filter(!is.na(yield), yield != 0)
    }
  )
}

crop_production <- map_dfr(
  crops,
  read_crop_production
)

crop_production <- crop_production %>%
  trim_xy()

saveRDS(
  crop_production,
  file.path(
    production_output_dir,
    "crop_production.rds"
  )
)