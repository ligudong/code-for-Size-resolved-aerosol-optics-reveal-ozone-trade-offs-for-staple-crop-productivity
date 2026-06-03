# ----------------------------
# 00_helper_functions.R
# Utility functions for raster checks, weighted R², formula generation,
# memory inspection, coordinate rounding, grid matching, data loading,
# and relative prediction with confidence intervals.
# ----------------------------

library(dplyr)
library(terra)
library(stringr)
library(arrow)
library(dtplyr)
library(readxl)
library(purrr)
library(tibble)
library(fixest)

# Check if object is a raster
is.Raster <- function(x) {
  class(x)[1] %in% c("RasterLayer","RasterBrick","RasterStack","SpatRaster")
}

# Weighted R²
get_wr2_ss <- function(y, y_pred, w) {
  ss_residual <- sum(w * (y - y_pred)^2)
  ss_total <- sum(w * (y - weighted.mean(y, w))^2)
  1 - ss_residual / ss_total
}

# Concatenate terms for formula
add_terms <- function(av) {
  reduce(av, function(x, y) str_c(x, " + ", y))
}

# Check memory usage
object_size <- function(aobject) {
  aobject %>%
    object.size() %>%
    print(unit = "auto")
}

# Set number of threads for qs
qn <- ifelse(Sys.info()[1]=="Windows",5,3)

# Round x,y coordinates to 4 decimals
trim_xy <- function(adata_xy) {
  adata_xy %>% mutate(across(c(x, y), ~round(.x,4)))
}

# Check if analysis data overlaps mask grid
check_join <- function(adata) {
  targeted <- rast("E:/data/outputs/masks/mask_extraction.tif") %>%
    sum(na.rm = T) %>%
    as.data.frame(xy = T) %>%
    distinct(x, y) %>%
    trim_xy()
  de <- adata %>% distinct(x, y)
  joined <- inner_join(targeted, de)
  antied <- anti_join(targeted, de)
  tibble(
    mask_grid = nrow(targeted),
    test_grid = nrow(de),
    inner_joined = nrow(joined),
    anti_joined = nrow(antied)
  )
}

# Default arrow for plots
plot_arrow <- arrow(length=unit(0.015,"npc"), ends="last", type="open")


get_data <- function(varss = NULL) {
  
  if (is.null(varss)) {
    varss <- c(
      "GOSIF_sum",
      "cloud",
      "Faod",
      "Caod",
      "AOT40",
      "W126",
      "O3",
      "O3_sum",
      "O3_peak",
      "VPD",
      "mean_surface",
      "mean_root",
      "maxtmp",
      "PM25",
      "PM10",
      "GPP",
      "irg_fraction",
      "RTSIF_sum",
      "CSIF_sum",
      "co2",
      str_c("bin", 1:42),
      str_c("step", 1:10),
      str_c("root_", 1:9),
      str_c("D", 1:13),
      str_c("d", 1:13)
    )
  }
  
  qread("E:/data/outputs/tidied.qs", nthreads = qn) %>%
    lazy_dt() %>%
    mutate(
      crop_parent = fifelse(str_detect(crop, "Rice"), "Rice", crop),
      x_y = str_c(crop, x_y),
      
      across(starts_with("GOSIF"), log),
      
      AOT40 = fcase(
        crop_parent == "Maize", AOT40_6,
        crop_parent == "Rice",  AOT40_7,
        crop_parent == "Wheat", AOT40_7 - AOT40_1
      ),
      
      W126 = fcase(
        crop_parent == "Maize", W126_6,
        crop_parent == "Rice",  W126_7,
        crop_parent == "Wheat", W126_7 - W126_1
      ),
      
      O3 = fcase(
        crop_parent == "Maize", O3_6,
        crop_parent == "Rice",  O3_7,
        crop_parent == "Wheat", (O3_7 * (MA - `GR&EM` + 1) - O3_1) / (MA - `GR&EM` + 0)
      ),
      
      O3_sum = fcase(
        crop_parent == "Maize", O3_sum_6,
        crop_parent == "Rice",  O3_sum_7,
        crop_parent == "Wheat", O3_sum_7 - O3_sum_1
      ),
      
      O3_peak = fcase(
        crop_parent == "Maize", O3_peak_6,
        crop_parent == "Rice",  O3_peak_7,
        crop_parent == "Wheat", pmax(O3_peak_1, O3_peak_2, O3_peak_3,
                                     O3_peak_4, O3_peak_5, O3_peak_6,
                                     O3_peak_7, na.rm = TRUE)
      )
    ) %>%
    drop_na(all_of(varss)) %>%
    group_by(crop, x, y) %>%
    filter(n() >= 10) %>%
    ungroup() %>%
    as_tibble() %>%
    inner_join(
      read_xlsx("E:/data/regions.xlsx") %>%
        mutate(region = str_c("R", region)),
      by = join_by(province, crop_parent)
    )
}

# Relative prediction with CI
relpred <- function(object, newdata, baseline=NULL, level=0.90){
  if(any(sapply(newdata, function(x) !is.numeric(x)))) stop("Only numeric variables supported")
  if(!is.null(baseline) & ncol(newdata)!=length(baseline)) stop("baseline length mismatch")
  newdata_filled <- fill_missing_vars(object,newdata)
  baseline_filled <- if(is.null(baseline)){ baseline_df <- newdata; baseline_df[] <- 0; fill_missing_vars(object,baseline_df) } else fill_missing_vars(object,as.data.frame(baseline))
  X <- as.matrix(newdata_filled - baseline_filled)
  B <- as.numeric(coef(object))
  df <- if(inherits(object,"fixest")) attributes(vcov(object, attr=T))$G else object$df.residual
  X <- X[,names(coef(object))]
  fit <- data.frame(fit=X%*%B)
  sig <- vcov(object)
  se <- apply(X,1,get_se,sig=sig)
  t_val <- qt((1-level)/2+level, df=df)
  fit$lwr <- fit$fit - t_val*se
  fit$upr <- fit$fit + t_val*se
  fit %>% as_tibble()
}

fill_missing_vars <- function(object,X){
  orig_vars <- all.vars(formula(object))[-1]
  fact_vars <- names(object$xlevels)
  for(f in fact_vars){
    if(!(f %in% names(X))) X[[f]] <- factor(object$xlevels[[f]][1], levels=object$xlevels[[f]])
    else X[[f]] <- factor(X[[f]], levels=object$xlevels[[f]])
  }
  for(v in orig_vars) if(!(v %in% names(X))) X[[v]] <- 0
  tt <- terms(object)
  Terms <- delete.response(tt)
  as.data.frame(model.matrix(Terms, data=X))
}

get_se <- function(r,sig){ sqrt(matrix(r,nrow=1) %*% sig %*% t(matrix(r,nrow=1))) }