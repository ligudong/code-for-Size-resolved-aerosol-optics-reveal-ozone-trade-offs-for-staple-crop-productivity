library(tidyverse)
library(terra)

# Define paths
data_file     <- "data/derived/panel_yucheng.rds"

rtsif_file    <- "data/derived/sif/rtsif.rds"
csif_file     <- "data/derived/sif/csif.rds"
gosif_file    <- "data/derived/sif/gosif.rds"

calendar_file <- "data/derived/calendar/tidied.rds"

output_file   <- "data/derived/panel_wheat_with_calendar.rds"

# Load Yucheng panel data
panel_data <- readRDS(data_file)

# Target location: Yucheng station
target_x <- 116 + 34 / 60 + 12.72 / 3600
target_y <- 36 + 49 / 60 + 44.4 / 3600

# Extract nearest SIF pixel time series
get_nearest_point <- function(df) {
  
  nearest_xy <- df %>%
    mutate(dist = sqrt((x - target_x)^2 + (y - target_y)^2)) %>%
    arrange(dist) %>%
    slice(1) %>%
    select(x, y)
  
  df %>%
    filter(x == nearest_xy$x, y == nearest_xy$y) %>%
    mutate(
      year = as.integer(year),
      month = as.integer(month)
    ) %>%
    select(-x, -y)
}

# Load and process RTSIF
RTSIF <- read_rds(rtsif_file) %>%
  mutate(RTSIF = rtsif * 0.0001)

RTSIF_point_sel <- get_nearest_point(RTSIF)

# Load and process CSIF
CSIF <- read_rds(csif_file) %>%
  mutate(CSIF = csif * 0.0001)

CSIF_point_sel <- get_nearest_point(CSIF)

# Load and process GOSIF
GOSIF <- read_rds(gosif_file) %>%
  mutate(GOSIF2 = GOSIF * 0.0001)

GOSIF_point_sel <- get_nearest_point(GOSIF)

# Load crop calendar and keep wheat growing-season calendar
cldr_wheat <- read_rds(calendar_file) %>%
  filter((MA - `GR&EM`) >= 2) %>%
  trim_xy() %>%
  filter(crop == "Wheat") %>%
  select(
    year,
    month,
    HE,
    MA,
    GR_EM = `GR&EM`
  )

# Merge three SIF products
panel_SIF <- RTSIF_point_sel %>%
  left_join(CSIF_point_sel, by = c("year", "month")) %>%
  left_join(GOSIF_point_sel, by = c("year", "month")) %>%
  select(-RTSIF, -CSIF)

# Merge SIF and crop calendar into Yucheng panel
panel_wheat <- panel_data %>%
  left_join(panel_SIF, by = c("年" = "year", "月" = "month")) %>%
  left_join(cldr_wheat, by = c("年" = "year", "月" = "month"))

# Save output
saveRDS(panel_wheat, output_file)

head(panel_wheat)
names(panel_wheat)