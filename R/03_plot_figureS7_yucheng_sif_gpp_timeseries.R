library(tidyverse)
library(lubridate)
library(hrbrthemes)

panel_file <- "data/derived/panel_wheat_with_calendar.rds"
figure_dir <- "figures/supplementary/figureS7"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

panel_wheat <- read_rds(panel_file) %>%
  distinct() %>%
  mutate(
    date = make_date(年, 月, 15),
    GPP = abs(GEE),
    GPP_scaled   = GPP / max(GPP, na.rm = TRUE),
    RTSIF_scaled = rtsif / max(rtsif, na.rm = TRUE),
    CSIF_scaled  = csif / max(csif, na.rm = TRUE),
    GOSIF_scaled = GOSIF / max(GOSIF, na.rm = TRUE),
    GR_EM_date = make_date(年, 1, 1) + days(GR_EM - 1),
    MA_date    = make_date(年, 1, 1) + days(MA - 1)
  )

p_s7a <- ggplot(panel_wheat, aes(x = date)) +
  geom_rect(
    data = panel_wheat %>% filter(month(date) %in% (month(GR_EM):month(MA))),
    aes(xmin = date - 15, xmax = date + 15, ymin = -Inf, ymax = Inf),
    fill = "grey90",
    alpha = 0.2
  ) +
  geom_line(aes(y = GPP_scaled, color = "GPP"), linewidth = 1.2) +
  geom_line(aes(y = RTSIF_scaled, color = "RTSIF"), linewidth = 1.2, linetype = "dotdash") +
  geom_line(aes(y = CSIF_scaled, color = "CSIF"), linewidth = 1.2) +
  geom_line(aes(y = GOSIF_scaled, color = "GOSIF"), linewidth = 1.2) +
  scale_color_manual(
    values = c(
      "GPP"   = "black",
      "RTSIF" = "#f38181",
      "CSIF"  = "#95e1d3",
      "GOSIF" = "#3d84a8"
    )
  ) +
  labs(
    x = "Date",
    y = "Scaled GPP / SIF",
    color = "Variable"
  ) +
  theme_ipsum_rc(base_size = 30, axis_title_size = 30) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30),
    panel.grid.major = element_line(color = "grey"),
    panel.grid.minor = element_blank()
  )

print(p_s7a)

ggsave(
  filename = file.path(figure_dir, "figureS7_yucheng_sif_gpp_timeseries.tif"),
  plot = p_s7a,
  width = 12,
  height = 5,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)