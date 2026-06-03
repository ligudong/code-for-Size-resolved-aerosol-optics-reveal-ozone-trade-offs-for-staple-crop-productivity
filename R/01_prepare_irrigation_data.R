source("R/00_load_packages.R")
source("R/00_helper_functions.R")

irrigation_raw_dir <- "data/raw/irrigation"
mask_dir <- "data/derived/masks"
irrigation_output_dir <- "data/derived/irrigation"

dir.create(irrigation_output_dir, recursive = TRUE, showWarnings = FALSE)

mask_path <- file.path(mask_dir, "mask_extraction.tif")
mask <- rast(mask_path)[[1]]

# 读取所有灌溉 raster 并重采样到 mask
a <- list.files(irrigation_raw_dir, full.names = TRUE) %>%
  pro_map(~ rast(.x) %>% resample(mask, method = "sum")) %>%
  rast()

# 读取基准 raster（年份 2000），空值设为 1
b <- rast(file.path(irrigation_raw_dir, "2000.tif"))
b[is.na(b)] <- 1
b <- resample(b, mask, method = "sum")

# 计算灌溉比例
a <- a / b
a[is.na(a)] <- 0

# 转为面板数据
irrigation_panel <- a %>%
  as.data.frame(xy = TRUE) %>%
  pivot_longer(
    -c(x, y),
    names_to = "year",
    names_transform = list(year = as.integer),
    values_to = "irg_fraction"
  ) %>%
  trim_xy()

check_join(irrigation_panel)

saveRDS(
  irrigation_panel,
  file.path(irrigation_output_dir, "tidied.rds")
)