# Load packages and helper functions
source("R/00_load_packages.R")
source("R/00_helper_functions.R")

# Directories
mask_dir <- "data/derived/masks"
o3_input_dir <- "data/raw/monthly_O3"
o3_output_dir <- "data/derived/ozone"
station_files <- c(
  "data/raw/station_data/China.rds",
  "data/raw/station_data/China_NS2021.rds"
)
dir.create(o3_output_dir, recursive = TRUE, showWarnings = FALSE)

# Load maize mask
mask <- rast(file.path(mask_dir, "mask_extraction.tif")) %>% sum(na.rm = TRUE)

# Load monthly O3 raster files and convert to dataframe
data <- list.files(file.path(o3_input_dir, "monthly_O3"), full.names = TRUE) %>%
  map(function(afile) {
    r <- rast(afile) %>%
      `names<-`(time(.))
    
    resample(r, mask, method = "near") %>%
      mask(mask) %>%
      as.data.frame(xy = TRUE) %>%
      pivot_longer(
        -c(x, y),
        names_to = "date",
        names_transform = list(date = ymd),
        values_to = "O3"
      )
  }, .progress = TRUE) %>%
  bind_rows() %>%
  filter(year(date) %in% 2005:2022) %>%
  trim_xy()

# Check grid consistency
check_join(data)

# Load station data and combine
data_station <- map_dfr(station_files, read_rds) %>%
  inner_join(by = c("station_id", "date")) %>%  # adjust join keys if necessary
  mutate(O3 = rollingo3 * 1e3 * 2)  # convert ppb to ug/m3

# Fit GAMs for W126 and AOT40 using parallel cluster
cl <- makeCluster(qn)

model_w126 <- bam(
  W126 ~ s(O3),
  data = data_station,
  family = tw(),
  chunk.size = 5000,
  cluster = cl
)

model_aot40 <- bam(
  AOT40 ~ s(O3),
  data = data_station,
  family = tw(),
  chunk.size = 5000,
  cluster = cl
)

# Predict W126 and AOT40 for raster data
data <- data %>%
  mutate(
    W126 = exp(predict(model_w126, ., cluster = cl)),
    AOT40 = exp(predict(model_aot40, ., cluster = cl)),
    year = year(date),
    month = month(date),
    .keep = "unused"
  )

# Save processed data and GAM models
saveRDS(data, file.path(o3_output_dir, "tidied.rds"))
saveRDS(list(model_aot40 = model_aot40, model_w126 = model_w126),
        file.path(o3_output_dir, "gam.rds"))

# Stop cluster
stopCluster(cl)