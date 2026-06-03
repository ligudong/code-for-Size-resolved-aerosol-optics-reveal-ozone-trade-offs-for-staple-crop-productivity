source("R/00_load_packages.R")
source("R/00_helper_functions.R")

province_plot_file <- "data/raw/shp/province.shp"
region_plot_file <- "data/raw/shp/province.shp"
ozone_rds_file <- "data/derived/ozone_summary_tables/impacts_ozone_summarised_base.rds"

figure_dir <- "figures/ozone"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

province_plot <- st_read(province_plot_file)
region_plot <- st_read(region_plot_file) %>%
  group_by(region) %>%
  summarise(geometry = st_union(geometry))

ozone <- read_rds(ozone_rds_file)

ozone_plot <- function(acrop, plot_legend=TRUE, subtitle=NULL) {
  p <- ozone %>%
    filter(crop_parent == acrop) %>%
    unnest_wider(ozone_results) %>%
    select(crop, year, crop_parent, province_level) %>%
    unnest() %>%
    group_by(crop_parent, 省, peak_level) %>%
    summarise(across(contains("%"), mean), .groups="drop") %>%
    filter(peak_level==100) %>%
    inner_join(province_plot, c("省"="name")) %>%
    ggplot() +
    geom_sf(aes(fill=`50%`*100, geometry=geometry), color=NA, size=0.2, show.legend=plot_legend) +
    geom_sf(data=region_plot, fill=NA, size=0.8) +
    scale_fill_gradientn(name="Percentage change in SIF", limits=c(-1,3),
                         colours=c("#1D39C4","#8BA2FF","#91D5FF","#BAE7FE","white","#FFCCD1","#FF929D","#E65E67","#B8292F","#800D00"),
                         na.value="grey90") +
    xlab(NULL) + ylab(NULL) +
    labs(tag=subtitle)
  return(p)
}

crop_list <- unique(ozone$crop_parent)
map_list <- map(crop_list, ~ ozone_plot(.x))
patch <- wrap_plots(map_list, ncol=2)
ggsave(file.path(figure_dir,"Ozone_Impact_Map.pdf"), patch, width=10, height=12)
ggsave(file.path(figure_dir,"Ozone_Impact_Map.tif"), patch, width=10, height=12, dpi=600)