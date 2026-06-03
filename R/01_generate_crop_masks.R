# Load helper functions and formula definitions
source("R/00_helper_functions.R")
source("R/00_formula_definitions.R")

library(terra)
library(purrr)
library(stringr)

# Define input and output directories
crop_mask_dir <- "data/raw/crop_masks"
mask_output_dir <- "data/derived/masks"
dir.create(mask_output_dir, showWarnings = FALSE, recursive = TRUE)

# Reference raster for extent
reference_raster <- rast(file.path(crop_mask_dir, "CHN_Maize_2019.tif"))
exts <- ext(reference_raster)

crops <- c("Maize", "Rice", "Wheat")

# Unified 0.5° mask for all major crops
all_files <- list.files(crop_mask_dir, full.names = TRUE) %>%
  map(rast) %>%
  map(~ extend(.x, exts)) %>%
  rast()

all_mask <- tapp(all_files, rep(crops, each = 20), "sum", na.rm = TRUE)
all_mask[all_mask > 0] <- 1
all_mask[all_mask == 0] <- NA
all_mask <- project(all_mask, "epsg:4326", method = "near")

degree_raster <- rast(resolution = 0.5) %>%
  crop(all_mask, snap = "out")

all_mask <- resample(all_mask, degree_raster, "sum")
writeRaster(all_mask, file.path(mask_output_dir, "mask_extraction.tif"), overwrite = TRUE)

# Annual 0.5° mask for each crop
for (acrop in crops) {
  crop_files <- list.files(crop_mask_dir, acrop, full.names = TRUE) %>%
    map(rast) %>%
    map(~ extend(.x, exts)) %>%
    rast()
  
  names(crop_files) <- 2000:2019
  crop_mask <- project(crop_files, "epsg:4326", method = "near")
  
  degree_raster <- rast(resolution = 0.5) %>%
    crop(crop_mask, snap = "out")
  
  crop_mask <- resample(crop_mask, degree_raster, "sum")
  crop_mask[crop_mask == 0] <- NA
  
  writeRaster(crop_mask, file.path(mask_output_dir, str_c("mask_", acrop, ".tif")), overwrite = TRUE)
}