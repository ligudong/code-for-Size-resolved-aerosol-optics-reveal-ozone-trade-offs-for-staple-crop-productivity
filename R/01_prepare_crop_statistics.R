source("R/00_helper_functions.R")
source("R/00_load_packages.R")

library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readxl)
library(readr)

crop_stats_dir <- "data/raw/crop_statistics"
crop_stats_output_dir <- "data/derived/crop_statistics"

area_dir <- file.path(crop_stats_dir, "area")
production_dir <- file.path(crop_stats_dir, "production")
population_file <- file.path(crop_stats_dir, "population.xlsx")

dir.create(crop_stats_output_dir, showWarnings = FALSE, recursive = TRUE)

read_crop_stat_files <- function(data_dir, value_name, crop_pattern) {
  tibble(
    file = list.files(data_dir, full.names = TRUE),
    data = map(file, ~ fread(.x, skip = 3, nrows = 31)),
    crop = map_chr(file, ~ fread(.x, skip = 1, nrows = 1) %>% pull(V1))
  ) %>%
    select(-file) %>%
    mutate(crop = str_remove_all(crop, crop_pattern)) %>%
    unnest(data) %>%
    rename(year = 时间) %>%
    pivot_longer(
      cols = -c(crop, year),
      names_to = "province",
      values_to = value_name
    ) %>%
    drop_na() %>%
    mutate(year = parse_number(year))
}

area_data <- read_crop_stat_files(
  data_dir = area_dir,
  value_name = "area",
  crop_pattern = "指标：|播种面积\\(千公顷\\)"
)

production_data <- read_crop_stat_files(
  data_dir = production_dir,
  value_name = "production",
  crop_pattern = "指标：|产量\\(万吨\\)"
)

population_data <- read_xlsx(population_file, range = "A4:AF26") %>%
  rename(year = 时间) %>%
  pivot_longer(
    cols = -year,
    names_to = "province",
    values_to = "population"
  ) %>%
  mutate(year = parse_number(year))

crop_statistics <- inner_join(
  area_data,
  production_data,
  by = c("crop", "province", "year")
) %>%
  inner_join(
    population_data,
    by = c("province", "year")
  ) %>%
  mutate(
    production = production * 1e4,
    crop = fcase(
      crop == "小麦", "Wheat",
      crop == "玉米", "Maize",
      crop == "稻谷", "Rice"
    )
  ) %>%
  drop_na() %>%
  group_by(crop, year) %>%
  summarise(
    yield = sum(production, na.rm = TRUE) / sum(area, na.rm = TRUE),
    area = sum(area, na.rm = TRUE),
    population = sum(population, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(
  crop_statistics,
  file.path(crop_stats_output_dir, "summary_by_crop_year.rds")
)