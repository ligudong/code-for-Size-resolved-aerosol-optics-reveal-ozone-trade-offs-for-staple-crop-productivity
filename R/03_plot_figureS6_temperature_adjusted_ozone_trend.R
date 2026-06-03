library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(terra)
library(openxlsx)
library(extrafont)
library(showtext)

showtext_auto()
loadfonts(device = "win")


#Load temperature data

data_tem <- readRDS("data/derived/temp/tmax.rds")  # Project-relative path

tem <- data_tem %>%
  rename(maize = lyr1, rice = lyr2, wheat = temp) %>%
  select(-files_)

data_long <- tem %>%
  pivot_longer(
    cols = c(maize, rice, wheat),
    names_to = "crop",
    values_to = "temp_value"
  )


#Extract O₃ raster values per year

tiff_folder <- "data/raw/ozone/O3_tif"  # Clean relative path
years <- 2001:2023
data_long$O3 <- NA

for(yr in years){
  cat("Processing year:", yr, "\n")
  
  rast_file <- file.path(tiff_folder, paste0("south", yr, ".tif"))
  o3_rast <- rast(rast_file)
  
  idx <- data_long$year == yr
  pts_vect <- vect(data_long[idx, c("x","y")], geom = c("x","y"), crs = crs(o3_rast))
  values <- terra::extract(o3_rast, pts_vect)[,2]
  data_long$O3[idx] <- values
}

# Compute temperature-adjusted O₃ residuals

data_long$O3_adj <- NA
crops <- unique(data_long$crop)

for(c in crops){
  idx <- data_long$crop == c & !is.na(data_long$O3) & !is.na(data_long$temp_value)
  temp_model <- lm(O3 ~ temp_value, data = data_long[idx, ], na.action = na.exclude)
  data_long$O3_adj[idx] <- resid(temp_model)
}

# Save adjusted ozone dataset
saveRDS(data_long, "data/derived/ozone/ozone_adjusted.rds")
data_long <- readRDS("data/derived/ozone/ozone_adjusted.rds")


#Compute annual mean residuals by crop

O3_trend_year <- data_long %>%
  group_by(crop, year) %>%
  summarise(mean_O3_adj = mean(O3_adj, na.rm = TRUE), .groups = "drop")

# Export annual residuals to Excel
write.xlsx(
  O3_trend_year,
  file = "data/derived/ozone/O3_adj_annual.xlsx",
  rowNames = FALSE
)


tiff("figures/supplementary/figureS6/O3_trend_annual.tif",
     width = 8, height = 6, units = "in", res = 300, compression = "lzw")

ggplot(O3_trend_year, aes(x = year, y = mean_O3_adj, color = crop)) +
  geom_line(size = 1) +
  geom_point() +
  labs(x = "Year", y = "O₃ residuals") +
  theme_minimal(base_family = "Arial", base_size = 50) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 16, face = "bold")
  )

dev.off()