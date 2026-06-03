# 00_load_packages.R

rm(list = ls())
gc()

library(magrittr, exclude = "set_names")

library(tidyverse)
library(data.table)
library(dtplyr)
library(collapse)
library(lubridate)
library(slider)
library(purrr)
library(purrrgress)

library(terra)
library(raster, exclude = "select")
library(sf)
library(exactextractr)
library(ncdf4)

library(readxl)
library(readr)
library(writexl)
library(openxlsx)
library(haven)
library(qs)

library(fixest)
library(broom)

library(tidymodels)
library(correlation)
library(Matrix)
library(matrixStats, exclude = "count")
library(splines)
library(mgcv)

library(ggplot2)
library(ggpmisc)
library(hrbrthemes)
library(cowplot)
library(ggh4x)
library(scales)
library(patchwork)
library(ggside)
library(ggdensity)
library(ggupset)
library(latex2exp)
library(ggrepel)
library(ggpointdensity)
library(viridis)
library(gghighlight)
library(biscale)
library(ggdist)
library(ggstar)
library(openair)
library(grid)

library(sysfonts)
library(showtext)
library(furrr)
library(lubridate)
library(mgcv)
library(parallel)



sf_use_s2(FALSE)

if (Sys.info()[1] == "Windows") {
  font_path <- "C:/Program Files/R/R-4.2.3/library/roboto-condensed/RobotoCondensed-Regular.ttf"
} else {
  font_path <- "C:/Program Files/R/R-4.2.3/library/roboto-condensed/RobotoCondensed-Regular.ttf"
}

if (file.exists(font_path)) {
  sysfonts::font_add("Roboto Condensed", font_path)
}

showtext::showtext_auto()

