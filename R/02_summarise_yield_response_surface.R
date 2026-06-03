source("R/00_load_packages.R")
source("R/00_helper_functions.R")

library(dplyr)
library(tidyr)
library(purrr)
library(terra)
library(qs)

# Directories
input_dir  <- "data/derived/impacts"     # 上一步计算的 fAOD/O3 情景结果
output_dir <- "data/derived/yield_response"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load yield panel data
yield_data <- readRDS("data/derived/combined_crop_production_2010_2020_data.rds")

# Load fAOD and O3 impacts (SIF response) from previous simulations
faod_impacts <- qread(file.path(input_dir, "impacts_Faod.qs"), nthreads = qn)
ozone_impacts <- qread(file.path(input_dir, "impacts_ozone_summarised_faod.rds"))

# Merge yield and pollution scenarios
yield_response <- map(c("Maize","Rice","Wheat"), function(acrop){
  
  faod_data <- faod_impacts %>%
    filter(crop_parent == acrop) %>%
    unnest(Faod_results) %>%
    select(x, y, year, faod_level, `50%`) %>%
    rename(faod_value = `50%`)
  
  ozone_data <- ozone_impacts %>%
    filter(crop_parent == acrop) %>%
    unnest(ozone_results) %>%
    select(x, y, year, peak_level, `50%`) %>%
    rename(o3_value = `50%`)
  
  # Cross join Faod and O3 levels
  combined <- expand.grid(
    x = unique(faod_data$x),
    y = unique(faod_data$y),
    faod_level = unique(faod_data$faod_level),
    peak_level = unique(ozone_data$peak_level)
  )
  
  combined <- combined %>%
    left_join(faod_data, by = c("x","y","faod_level")) %>%
    left_join(ozone_data, by = c("x","y","peak_level")) %>%
    mutate(combined_sif = faod_value + o3_value) %>%
    select(x,y,faod_level,peak_level,combined_sif)
  
  combined
})

names(yield_response) <- c("Maize","Rice","Wheat")

# Save combined yield response
qsave(yield_response, file.path(output_dir, "yield_response_surface.rds"), nthread = qn)