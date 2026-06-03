# Load packages and helper functions
source("R/00_load_packages.R")
source("R/00_helper_functions.R")

# Directories
calendar_file <- "data/derived/calendar/tidied.rds"
fraction_dir  <- "data/derived/masks"
faod_file     <- "data/derived/aerosol/Faod.rds"
ozone_file    <- "data/derived/ozone/tidied.rds"
faod_cft_file <- "data/derived/faod_scenario/faod_cft_scenario.rds"
output_file   <- "data/derived/impacts/impacts_Faod.qs"

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Load crop calendar and filter valid growth periods
cldr <- read_rds(calendar_file) %>%
  filter((MA - `GR&EM`) >= 2) %>%
  trim_xy()

# Load crop fraction masks
fraction <- tibble(
  crop = c("Maize", "Rice(LR)", "Rice(SR&ER)", "Wheat"),
  data = map(c("Maize", "Rice", "Rice", "Wheat"), function(acrop) {
    afile <- str_c(fraction_dir, "/mask_", acrop, ".tif")
    rast(afile) %>%
      as.data.frame(xy = TRUE, na.rm = FALSE) %>%
      pivot_longer(
        -c(x, y),
        names_to = "year",
        names_transform = list(year = as.integer),
        values_drop_na = TRUE,
        values_to = "fraction"
      )
  })
) %>%
  unnest() %>%
  trim_xy()

# Load bootstrap regression results
f1 <- read_rds(file.path("data/derived/bootstrap", "boots_f3.rds")) %>%
  mutate(crop = crop_parent) %>%
  bind_rows((.) %>% filter(crop == "Rice") %>% mutate(crop = "Rice(LR)")) %>%
  mutate(crop = fifelse(crop == "Rice", "Rice(SR&ER)", crop))

# Load aerosol and ozone data
Faod <- read_rds(faod_file)
ozone <- read_rds(ozone_file)

# Remove rows with NA O3
rows_to_remove <- ozone %>%
  filter(is.na(O3)) %>%
  select(x, y, year) %>% distinct()

ozone <- ozone %>% anti_join(rows_to_remove, by = c("x", "y", "year"))

# Merge Faod and ozone data
Faod <- merge(
  ozone[, c("x", "y", "year", "month", "AOT40")],
  Faod[, c("x", "y", "year", "month", "Faod")],
  by = c("x", "y", "year", "month"), all = TRUE
) %>% drop_na()

# Aggregate monthly mean per crop
Faod <- cldr %>%
  lazy_dt() %>%
  inner_join(Faod) %>%
  group_by(crop, x, y, year) %>%
  summarise(Faod = mean(Faod),
            AOT40 = mean(AOT40), .groups = "drop") %>%
  inner_join(fraction) %>%
  as_tibble() %>%
  nest(fdata = -c(crop, year)) %>%
  arrange(crop, year)

# Load Faod counterfactual scenarios
faod_ctr <- read_rds(faod_cft_file)
faod_ctr <- faod_ctr %>%
  mutate(faod_data = map(faod_data, function(adata) {
    cldr %>%
      lazy_dt() %>%
      inner_join(adata, by = c("x", "y", "month")) %>%
      group_by(crop, year, x, y) %>%
      summarise(Faod_cft = mean(Faod_cft), .groups = "drop") %>%
      as_tibble()
  }, .progress = TRUE)) %>%
  unnest() %>%
  nest(faod_data = -c(crop, year))

# Merge observed Faod, coefficients, and counterfactuals
impacts <- reduce(list(Faod, f1, faod_ctr), inner_join)

# Parallel computation for relative Faod impacts
plan(multisession, workers = 6)
impacts <- impacts %>%
  mutate(Faod_results = future_pmap(
    list(fdata, coefs, faod_data),
    function(adata, coefs, faod_data) {
      
      coefs <- coefs %>%
        select(id, term, estimate) %>%
        pivot_wider(names_from = term, values_from = estimate) %>%
        select(contains("Faod"))
      
      cal_data <- adata %>%
        inner_join(faod_data, by = c("x", "y")) %>%
        mutate(Faod_cft = fifelse(Faod_cft > Faod, Faod, Faod_cft))
      
      rel_X <- cal_data %>%
        mutate(Faod1 = Faod_cft - Faod,
               Faod2 = Faod_cft^2 - Faod^2,
               OFaod1 = Faod1 * AOT40,
               OFaod2 = Faod2 * AOT40) %>%
        select(Faod1, Faod2, OFaod1, OFaod2)
      
      rel_results <- tcrossprod(as.matrix(rel_X), as.matrix(coefs)) %>%
        expm1() %>%
        rowQuantiles(probs = c(0.05, 0.5, 0.95)) %>%
        as_tibble() %>%
        bind_cols(cal_data %>% select(x, y, fraction, faod_level))
      
      return(rel_results)
    }, .progress = TRUE
  )) %>%
  select(-c(fdata, coefs, faod_data))

plan(sequential)

# Example plot for debugging
impacts$Faod_results[[1]] %>%
  filter(faod_level == 0.40) %>%
  select(x, y, `50%`) %>%
  rast(type = "xyz") %>%
  plot()

# Save final relative impact results
qsave(impacts, output_file, nthread = qn)