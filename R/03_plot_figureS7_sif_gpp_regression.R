library(tidyverse)
library(fixest)
library(broom)
library(openxlsx)

panel_file <- "data/derived/panel_wheat_with_calendar.rds"
figure_dir <- "figures/supplementary/figureS7"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

panel_wheat <- read_rds(panel_file) %>%
  distinct()

panel_wheat_season <- panel_wheat %>%
  mutate(
    GPP = -GEE
  ) %>%
  drop_na(
    GPP,
    GOSIF,
    csif,
    rtsif,
    `月平均近地面空气温度`,
    `月平均近地面空气湿度`,
    `月平均光合有效辐射`,
    `风向`,
    `大气压`,
    `月降水量`
  ) %>%
  mutate(
    GPP_scaled   = as.numeric(scale(GPP)),
    GOSIF_scaled = as.numeric(scale(GOSIF)),
    CSIF_scaled  = as.numeric(scale(csif)),
    RTSIF_scaled = as.numeric(scale(rtsif)),
    temp_scaled  = as.numeric(scale(`月平均近地面空气温度`)),
    hum_scaled   = as.numeric(scale(`月平均近地面空气湿度`)),
    par_scaled   = as.numeric(scale(`月平均光合有效辐射`)),
    wind_scaled  = as.numeric(scale(`风向`)),
    pres_scaled  = as.numeric(scale(`大气压`)),
    prec_scaled  = as.numeric(scale(`月降水量`))
  )

formulas <- list(
  RTSIF = GPP_scaled ~ RTSIF_scaled +
    temp_scaled + hum_scaled + par_scaled +
    wind_scaled + pres_scaled + prec_scaled,
  
  GOSIF = GPP_scaled ~ GOSIF_scaled +
    temp_scaled + hum_scaled + par_scaled +
    wind_scaled + pres_scaled + prec_scaled,
  
  CSIF = GPP_scaled ~ CSIF_scaled +
    temp_scaled + hum_scaled + par_scaled +
    wind_scaled + pres_scaled + prec_scaled
)

models <- lapply(formulas, function(fml) {
  feols(fml, data = panel_wheat_season)
})

result_table <- bind_rows(
  lapply(names(models), function(m) {
    tidy(models[[m]]) %>%
      mutate(model = m)
  })
) %>%
  mutate(
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    )
  )

model_stats <- bind_rows(
  lapply(names(models), function(m) {
    tibble(
      model = m,
      r2 = as.numeric(fitstat(models[[m]], "r2")),
      adj_r2 = as.numeric(fitstat(models[[m]], "ar2")),
      n = nobs(models[[m]])
    )
  })
)

result_export <- result_table %>%
  select(
    model,
    term,
    estimate,
    std.error,
    statistic,
    stars
  ) %>%
  mutate(
    term = gsub("_scaled", "", term)
  ) %>%
  rename(
    Model = model,
    Variable = term,
    Estimate = estimate,
    Std_Error = std.error,
    t_value = statistic,
    Significance = stars
  )

write.xlsx(
  list(
    "Regression_results" = result_export,
    "Model_statistics" = model_stats
  ),
  file = file.path(figure_dir, "figureS7_sif_gpp_regression_results.xlsx"),
  rowNames = FALSE
)

print(result_export)
print(model_stats)