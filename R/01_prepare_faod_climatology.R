source("R/00_load_packages.R")
source("R/00_helper_functions.R")

# Directories
aerosol_dir <- "data/derived/aerosol"
output_dir  <- "data/derived/faod_scenario"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load Faod panel data
faod <- read_rds(file.path(aerosol_dir, "Faod.rds"))
summary(faod)

# Compute monthly climatology
month_climatology <- faod %>%
  fgroup_by(x, y, month) %>%                       # Group by location and month
  fsummarise(Faod_month = fmean(Faod)) %>%        # Compute mean Faod per month
  fgroup_by(x, y) %>%                             # Regroup by location
  fmutate(Faod_month_percent = Faod_month / fmean(Faod_month), .keep = "unused") %>% # Relative monthly percentage
  fungroup()

# Generate Faod counterfactual for different pollution levels
faod_cft_all <- tibble(
  faod_level = seq(0, 2.75, by = 0.01),  # Sequence of pollution levels
  faod_data = map(faod_level, function(anum) {
    month_climatology %>%
      distinct(x, y) %>%                 # Keep unique grid points
      mutate(Faod = pmax(anum, 0)) %>%   # Ensure Faod >= 0
      inner_join(month_climatology) %>%  # Join with monthly climatology
      mutate(Faod_cft = Faod_month_percent * Faod, .keep = "unused") # Compute adjusted Faod
  })
)

# Save counterfactual Faod scenarios
saveRDS(faod_cft_all, file.path(output_dir, "faod_cft_scenario.rds"))