source("R/00_load_packages.R")
source("R/00_helper_functions.R")

input_file <- "data/derived/yield_response/yield_response_scenarios.rds"
output_dir <- "figures/p5"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

yield_response <- readRDS(input_file)


p5c <- yield_response %>%
  filter(faod_level == 0.40, peak_level == 100) %>%
  group_by(year) %>%
  summarise(across(contains("kcal"), sum)) %>%
  ggplot(aes(x = year)) +
  geom_hline(aes(yintercept = kcal_percapita_perday), linetype="longdash", color="grey30") +
  geom_hline(aes(yintercept = kcal_percapita_perday), data=. %>% filter(year==2015), linetype="dashed", color="grey50") +
  geom_col(aes(y = kcal_percapita_perday_our), fill="#F9A828", width=0.5, color="black") +
  geom_errorbar(aes(ymax = kcal_percapita_perday_our_up, ymin = kcal_percapita_perday_our_lw), width=0.25) +
  geom_col(aes(y = kcal_percapita_perday), fill="#009d73", width=0.5, color="black") +
  scale_y_continuous(name="Calorie intake (kCal capita⁻¹ day⁻¹)") +
  scale_x_continuous(name=NULL, breaks=c(2005,2010,2015,2019)) +
  coord_cartesian(ylim=c(850,NA)) +
  theme_minimal(base_size=16)

ggsave(file.path(output_dir,"figure5c.tif"), plot=p5c, width=8, height=4, dpi=300)



# Example data (replace with real processed scenario if available)
data_5b <- data.frame(
  Scenario = rep("Rice",3),
  O3_Level = c("O₃ = 60","O₃ = 80","O₃ = 100"),
  value = c(16.32,30.17,39.04)
)
data_5b$O3_Level <- factor(data_5b$O3_Level, levels=c("O₃ = 100","O₃ = 80","O₃ = 60"))

p5b <- ggplot(data_5b, aes(x=Scenario, y=value)) +
  geom_bar(stat="identity", fill="#B7EB8F", width=0.6, color="black") +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  labs(y="Ozone-induced crop production change (Tg/yr)",
       title="Ozone Impact on Crop Yield by Scenario") +
  theme_minimal(base_size=14) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=14, face="bold", color="black"),
        axis.title.y=element_text(size=16, face="bold"),
        plot.title=element_text(size=18, hjust=0.5, face="bold"),
        legend.position="none")

ggsave(file.path(output_dir,"figure5b.tif"), plot=p5b, width=10, height=6, dpi=300)