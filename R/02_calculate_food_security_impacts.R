source("R/00_load_packages.R")
source("R/00_helper_functions.R")

# Directories
input_dir <- "data/derived"
output_dir <- "data/derived/yield_response"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load area, production, population
Area <- read_xlsx(file.path("data/raw", "area.xlsx"))
Production <- read_xlsx(file.path("data/raw", "production.xlsx"))
Pop <- read_xlsx(file.path("data/raw", "pop.xlsx"))

# Combine data
stats <- Area %>%
  inner_join(Production, by = c("year","crop","province")) %>%
  inner_join(Pop, by = c("year","province")) %>%
  mutate(production = production * 1e4) %>%  # kg
  drop_na() %>%
  group_by(crop, year) %>%
  summarise(
    yield = sum(production)/sum(area),
    area = sum(area),
    population = sum(population),
    .groups = "drop"
  )

# Load SIF response scenario data
Faod_full <- qread(file.path(input_dir,"impacts_Faod.qs"), nthreads = qn)
ozone_full <- qread(file.path(input_dir,"impacts_ozone_summarised_faod.rds"))

# Join datasets and calculate per-capita calorie intake
yield_response <- stats %>%
  inner_join(Faod_full, by = c("crop" = "crop_parent", "year")) %>%
  inner_join(ozone_full, by = c("crop" = "crop_parent", "year")) %>%
  mutate(
    n = fcase(crop=="Wheat", 0.78, crop=="Rice",1, crop=="Maize",0.79),
    w = fcase(crop=="Wheat", 0.2, crop=="Rice",0.1, crop=="Maize",0.7),
    E = fcase(crop=="Wheat", 3391.67, crop=="Rice",3882.05, crop=="Maize",3622.95),
    kcal_percapita_perday_our = 0.44 * area * (yield + `aer50%`*yield + `ozo50%`*yield) * 1e3 * n * (1-w) * E / (population*1e4)/365,
    kcal_percapita_perday_our_up = 0.44 * area * (yield + `aer95%`*yield + `ozo95%`*yield) * 1e3 * n * (1-w) * E / (population*1e4)/365,
    kcal_percapita_perday_our_lw = 0.44 * area * (yield + `aer5%`*yield + `ozo5%`*yield) * 1e3 * n * (1-w) * E / (population*1e4)/365
  )

# Save processed data for plotting
saveRDS(yield_response, file.path(output_dir,"yield_response_scenarios.rds"))