source("R/00_load_packages.R")
source("R/00_helper_functions.R")

vpd_raw_dir <- "data/raw/era5_climate"
mask_dir <- "data/derived/masks"
vpd_output_dir <- "data/derived/vpd"

dir.create(
  vpd_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

mask <- rast(
  file.path(mask_dir, "mask_extraction.tif")
) %>%
  sum(na.rm = TRUE)

VPD <- tibble(
  year = 2005:2019,
  data = map(
    year,
    function(ayear) {
      
      pressure <- rast(
        file.path(
          vpd_raw_dir,
          str_c(
            "e5.moda.an.sfc.128_151_msl.ll025sc.",
            ayear,
            "010100_",
            ayear,
            "120100.nc"
          )
        )
      ) / 100
      
      temperature <- rast(
        file.path(
          vpd_raw_dir,
          str_c(
            "e5.moda.an.sfc.128_167_2t.ll025sc.",
            ayear,
            "010100_",
            ayear,
            "120100.nc"
          )
        )
      ) - 273.15
      
      dewpoint <- rast(
        file.path(
          vpd_raw_dir,
          str_c(
            "e5.moda.an.sfc.128_168_2d.ll025sc.",
            ayear,
            "010100_",
            ayear,
            "120100.nc"
          )
        )
      ) - 273.15
      
      pressure <- pressure %>%
        `crs<-`("+proj=longlat +datum=WGS84") %>%
        project(mask) %>%
        mask(mask)
      
      temperature <- temperature %>%
        `crs<-`("+proj=longlat +datum=WGS84") %>%
        project(mask) %>%
        mask(mask)
      
      dewpoint <- dewpoint %>%
        `crs<-`("+proj=longlat +datum=WGS84") %>%
        project(mask) %>%
        mask(mask)
      
      fw <- 1 + 7e-4 + 3.46 * 1e-6 * pressure
      
      svp <- 6.112 * fw * exp(
        17.67 * temperature / (temperature + 243.5)
      )
      
      avp <- 6.112 * fw * exp(
        17.67 * dewpoint / (dewpoint + 243.5)
      )
      
      vpd <- svp - avp
      names(vpd) <- 1:12
      
      vpd %>%
        as.data.frame(xy = TRUE)
    },
    .progress = TRUE
  )
)

VPD <- VPD %>%
  unnest(data) %>%
  pivot_longer(
    cols = -c(x, y, year),
    names_to = "month",
    names_transform = list(month = as.integer),
    values_to = "VPD"
  ) %>%
  trim_xy()

check_join(VPD)

saveRDS(
  VPD,
  file.path(vpd_output_dir, "tidied.rds")
)