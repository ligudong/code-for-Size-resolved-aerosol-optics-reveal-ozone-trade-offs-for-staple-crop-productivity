source("R/00_load_packages.R")
source("R/00_helper_functions.R")

csif_dir <- "data/raw/csif"
mask_dir <- "data/derived/masks"
sif_output_dir <- "data/derived/sif"

dir.create(
  sif_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

mask <- rast(
  file.path(mask_dir, "mask_extraction.tif")
)[[1]]

nc_files <- list.files(
  csif_dir,
  pattern = "\\.nc$",
  full.names = TRUE
)

CSIF <- map_dfr(
  nc_files,
  function(nc_file) {
    
    year_current <- str_extract(
      basename(nc_file),
      "\\d{4}"
    ) %>%
      as.integer()
    
    nc_data <- rast(nc_file)
    
    map_dfr(
      seq_len(nlyr(nc_data)),
      function(i) {
        
        r <- nc_data[[i]]
        
        resample(
          r,
          mask,
          method = "near"
        ) %>%
          mask(mask) %>%
          as.data.frame(
            xy = TRUE,
            na.rm = FALSE
          ) %>%
          rename(
            CSIF = 3
          ) %>%
          mutate(
            year = year_current,
            month = i
          ) %>%
          select(
            x,
            y,
            year,
            month,
            CSIF
          )
      }
    )
  }
)

CSIF <- CSIF %>%
  trim_xy()

check_join(CSIF)

saveRDS(
  CSIF,
  file.path(
    sif_output_dir,
    "CSIF.rds"
  )
)