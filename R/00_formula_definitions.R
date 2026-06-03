
fml_base <- str_c("GOSIF_sum ~
            Caod + Caod ^ 2 +
            Faod + Faod ^ 2 +
            cloud + cloud ^ 2 + AOT40 + Faod*AOT40 + (Faod^ 2)*AOT40 + co2+
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()


fml_base2 <- str_c("CSIF_sum ~
            Caod + Caod ^ 2 +
            Faod + Faod ^ 2 +
            cloud + cloud ^ 2 + AOT40 + Faod*AOT40 + (Faod^ 2)*AOT40 + co2+
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()

fml_base3 <- str_c("RTSIF_sum ~
            Caod + Caod ^ 2 +
            Faod + Faod ^ 2 +
            cloud + cloud ^ 2 + AOT40 + Faod*AOT40 + (Faod^ 2)*AOT40 + co2+
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()

fml_base4 <- str_c("GOSIF_sum ~
            Caod + Caod ^ 2 +
            Faod + Faod ^ 2 +
            cloud + cloud ^ 2 + AOT40 + Faod*AOT40 + (Faod^ 2)*AOT40 + 
             + Caod*AOT40 + (Caod^ 2)*AOT40 + co2 +
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()


fml_base_no_interaction <- str_c("GOSIF_sum ~
            Caod + Caod ^ 2 +
            Faod + Faod ^ 2 +
            cloud + cloud ^ 2 + AOT40 + co2 +
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()

fml_base_no_interaction2 <- str_c("CSIF_sum ~
            Caod + Caod ^ 2 +
            Faod + Faod ^ 2 +
            cloud + cloud ^ 2 + AOT40 + co2 +
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()

fml_base_no_interaction3 <- str_c("RTSIF_sum ~
            Caod + Caod ^ 2 +
            Faod + Faod ^ 2 +
            cloud + cloud ^ 2 + AOT40 + co2 +
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()

fml_base_no_interaction4 <- str_c("RTSIF_sum ~
            Caod + Caod ^ 2 +
            Faod + Faod ^ 2 +
            cloud + cloud ^ 2 + AOT40 + co2 + + Faod*AOT40 + (Faod^ 2)*AOT40 +
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()

fml_base_no_interaction4 <- str_c("RTSIF_sum ~
            Caod + Caod ^ 2 +
            Faod + Faod ^ 2 +
            cloud + cloud ^ 2 + AOT40 + co2 +
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()


fml_gosif_aot40_only <- str_c("GOSIF_sum ~ AOT40 | x_y[year]") %>%
  as.formula()

fml_main_linear_interaction <- str_c("GOSIF_sum ~
            Caod + Faod + cloud  + AOT40 + Faod*AOT40 +
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()

#######irrigation heterogeneity interaction model############
fml_irrigation_heterogeneity_gosif <- str_c(
  "GOSIF_sum ~ ",
  c(
    str_c(c(
      "Faod", "Faod ^ 2","Caod", "Caod ^ 2", "AOT40", "cloud", "cloud ^ 2",
      "Faod*AOT40", "(Faod^ 2)*AOT40", "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9)
    ), " * irg_fraction")
  ) %>% str_flatten(collapse = "+"), "| x_y[year]"
) %>%
  as.formula()

fml_irrigation_heterogeneity_csif <- str_c(
  "CSIF_sum ~ ",
  c(
    str_c(c(
      "Faod", "Faod ^ 2","Caod", "Caod ^ 2", "AOT40", "cloud", "cloud ^ 2",
      "Faod*AOT40", "(Faod^ 2)*AOT40", "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9)
    ), " * irg_fraction")
  ) %>% str_flatten(collapse = "+"), "| x_y[year]"
) %>%
  as.formula()

fml_irrigation_heterogeneity_rtsif <- str_c(
  "RTSIF_sum ~ ",
  c(
    str_c(c(
      "Faod", "Faod ^ 2","Caod", "Caod ^ 2", "AOT40", "cloud", "cloud ^ 2",
      "Faod*AOT40", "(Faod^ 2)*AOT40", "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9)
    ), " * irg_fraction")
  ) %>% str_flatten(collapse = "+"), "| x_y[year]"
) %>%
  as.formula()

#######irrigation heterogeneity model without Faod × AOT40 interaction######
fml_irrigation_heterogeneity_no_interaction_gosif <- str_c(
  "GOSIF_sum ~ ",
  c(
    str_c(c(
      "Faod", "Faod ^ 2",
      "Caod", "Caod ^ 2",
      "AOT40",
      "cloud", "cloud ^ 2",
      "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9)
    ), " * irg_fraction")
  ) %>% str_flatten(collapse = "+"),
  "| x_y[year]"
) %>%
  as.formula()

fml_irrigation_heterogeneity_no_interaction_csif <- str_c(
  "CSIF_sum ~ ",
  c(
    str_c(c(
      "Faod", "Faod ^ 2",
      "Caod", "Caod ^ 2",
      "AOT40",
      "cloud", "cloud ^ 2",
      "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9)
    ), " * irg_fraction")
  ) %>% str_flatten(collapse = "+"),
  "| x_y[year]"
) %>%
  as.formula()

fml_irrigation_heterogeneity_no_interaction_rtsif <- str_c(
  "RTSIF_sum ~ ",
  c(
    str_c(c(
      "Faod", "Faod ^ 2",
      "Caod", "Caod ^ 2",
      "AOT40",
      "cloud", "cloud ^ 2",
      "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9)
    ), " * irg_fraction")
  ) %>% str_flatten(collapse = "+"),
  "| x_y[year]"
) %>%
  as.formula()


###### PM10/PM25 base model #####
fml_pm_interaction_gosif <- str_c(
  "GOSIF_sum ~
   PM10 + I(PM10^2) +
   PM25 + I(PM25^2) +
   cloud + I(cloud^2) +
   AOT40 + PM25*AOT40 + I(PM25^2)*AOT40 + co2 +
   ",
  add_terms(c(
    str_c("bin", 1:42),
    str_c("root_", 1:9),
    str_c("d", 1:13)
  )),
  "| x_y[year]"
) %>%
  as.formula()

fml_pm_base_csif <- str_c(
  "CSIF_sum ~
   PM10 + I(PM10^2) +
   PM25 + I(PM25^2) +
   cloud + I(cloud^2) +
   AOT40 + PM25*AOT40 + I(PM25^2)*AOT40 + co2 +
   ",
  add_terms(c(
    str_c('bin', 1:42),
    str_c('root_', 1:9),
    str_c('d', 1:13)
  )),
  "| x_y[year]"
) %>%
  as.formula()

fml_pm_base_rtsif <- str_c(
  "RTSIF_sum ~
   PM10 + I(PM10^2) +
   PM25 + I(PM25^2) +
   cloud + I(cloud^2) +
   AOT40 + PM25*AOT40 + I(PM25^2)*AOT40 + co2 +
   ",
  add_terms(c(
    str_c('bin', 1:42),
    str_c('root_', 1:9),
    str_c('d', 1:13)
  )),
  "| x_y[year]"
) %>%
  as.formula()



# Faod-only specification
fml_faod <- as.formula(
  "GOSIF_sum ~ Faod + Faod ^ 2 | x_y[year]"
)


# AOT40-only specification
fml_aot40 <- str_c(
  "GOSIF_sum ~
   AOT40 | x_y[year]"
) %>%
  as.formula()


# Additive specification with Faod and AOT40
fml_af <- str_c(
  "GOSIF_sum ~
   AOT40 + Faod + Faod ^ 2 |
   x_y[year]"
) %>%
  as.formula()


# Interaction specification between Faod and AOT40
fml_faod_aot40 <- str_c(
  "GOSIF_sum ~
   Faod + Faod ^ 2 +
   AOT40 +
   Faod * AOT40 + co2 +
   (Faod ^ 2) * AOT40 |
   x_y[year]"
) %>%
  as.formula()

#######W126/Ozone######
fml_base_O3 <- str_c(
  "GOSIF_sum ~
   Caod + I(Caod^2) +
   Faod + I(Faod^2) +
   cloud + I(cloud^2) +
   O3 + Faod*O3 + I(Faod^2)*O3 +co2 +
   ",
  add_terms(c(
    str_c("bin", 1:42),
    str_c("root_", 1:9),
    str_c("d", 1:13)
  )),
  "| x_y[year]"
) %>%
  as.formula()

fml_o3_no_interaction <- str_c("GOSIF_sum ~
            Caod + Caod ^ 2 +
            Faod + Faod ^ 2 +
            cloud + cloud ^ 2 + O3 + co2 +
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()

fml_base_W126 <- str_c(
  "GOSIF_sum ~
   Caod + I(Caod^2) +
   Faod + I(Faod^2) +
   cloud + I(cloud^2) +
   W126 + Faod*W126 + I(Faod^2)*W126 + co2 +
   ",
  add_terms(c(
    str_c("bin", 1:42),
    str_c("root_", 1:9),
    str_c("d", 1:13)
  )),
  "| x_y[year]"
) %>%
  as.formula()


fml_W126_no_interaction <- str_c(
  "GOSIF_sum ~
   Caod + I(Caod^2) +
   Faod + I(Faod^2) +
   cloud + I(cloud^2) +
   W126 + co2 +
   ",
  add_terms(c(
    str_c("bin", 1:42),
    str_c("root_", 1:9),
    str_c("d", 1:13)
  )),
  "| x_y[year]"
) %>%
  as.formula()

# Region-specific heterogeneous model across multiple regions
fml_inter_region <- str_c(
  "GOSIF_sum ~ ",
  c(
    str_c("R", 1:6, "_Caod"),
    str_c("R", 1:6, "_Caod2"),
    str_c("R", 1:6, "_Faod"),
    str_c("R", 1:6, "_Faod2"),
    str_c("R", 1:4, "_AOT40"),
    "cloud", "I(cloud^2)","co2",
    str_c("bin", 1:42),
    str_c("root_", 1:9),
    str_c("d", 1:13)
  ) %>%
    str_flatten(collapse = "+"),
  "| x_y[year]"
) %>%
  as.formula()

fml_no_Faod <- str_c("GOSIF_sum ~
            Caod + Caod ^ 2 +
            cloud + cloud ^ 2 + AOT40 + 
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()

fml_no_O3 <- str_c("GOSIF_sum ~
            Caod + Caod ^ 2 +
            Faod + Faod ^ 2 +
            cloud + cloud ^ 2 + 
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()

fml_no_Caod <- str_c("GOSIF_sum ~
            Faod + Faod ^ 2 +
            cloud + cloud ^ 2 + AOT40 + Faod*AOT40 +  (Faod^ 2)*AOT40 +
            ", add_terms(c(
              str_c("bin", 1:42),
              str_c("root_", 1:9),
              str_c("d", 1:13)
            )), "|
            x_y[year]") %>%
  as.formula()

# Irrigation heterogeneity model for GOSIF with CO2 and dummy controls
fml_irrigation_heterogeneity_gosif <- str_c(
  "GOSIF_sum ~ ",
  c(
    str_c(c(
      "Faod", "Faod ^ 2",
      "Caod", "Caod ^ 2",
      "AOT40", "cloud", "cloud ^ 2",
      "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9),
      str_c("d", 1:13)
    ), " * irg_fraction")
  ) %>% str_flatten(collapse = "+"),
  "| x_y[year]"
) %>%
  as.formula()

# Year-specific heterogeneous model for GOSIF with CO2 and dummy controls
fml_year_heterogeneity_gosif <- str_c(
  "GOSIF_sum ~ ",
  c(
    str_c("new_", c(
      "Faod", "Faod ^ 2",
      "Caod", "Caod ^ 2",
      "AOT40",
      "cloud", "cloud ^ 2",
      "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9),
      str_c("d", 1:13)
    )),
    str_c("mid_", c(
      "Faod", "Faod ^ 2",
      "Caod", "Caod ^ 2",
      "AOT40",
      "cloud", "cloud ^ 2",
      "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9),
      str_c("d", 1:13)
    )),
    str_c("old_", c(
      "Faod", "Faod ^ 2",
      "Caod", "Caod ^ 2",
      "AOT40",
      "cloud", "cloud ^ 2",
      "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9),
      str_c("d", 1:13)
    ))
  ) %>%
    str_flatten(collapse = "+"),
  "| x_y[year]"
) %>%
  as.formula()

# Irrigation heterogeneity model for GOSIF with CO2, dummy controls, and VPD interaction
fml_irrigation_heterogeneity_gosif <- str_c(
  "GOSIF_sum ~ ",
  c(
    str_c(c(
      "Faod", "Faod ^ 2",
      "Faod:AOT40", "I(Faod ^ 2):AOT40",
      "Caod", "Caod ^ 2",
      "AOT40", "cloud", "cloud ^ 2",
      "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9),
      str_c("d", 1:13)
    ), " * irg_fraction * VPD")
  ) %>% str_flatten(collapse = "+"),
  "| x_y[year]"
) %>%
  as.formula()

fml_irrigation_heterogeneity_gosif_no_interaction <- str_c(
  "GOSIF_sum ~ ",
  c(
    str_c(c(
      "Faod", "Faod ^ 2",
      "Caod", "Caod ^ 2",
      "AOT40", "cloud", "cloud ^ 2",
      "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9),
      str_c("d", 1:13)
    ), " * irg_fraction * VPD")
  ) %>% str_flatten(collapse = "+"),
  "| x_y[year]"
) %>%
  as.formula()

fml_irrigation_heterogeneity_csif <- str_c(
  "CSIF_sum ~ ",
  c(
    str_c(c(
      "Faod", "Faod ^ 2",
      "Faod:AOT40", "I(Faod ^ 2):AOT40",
      "Caod", "Caod ^ 2",
      "AOT40", "cloud", "cloud ^ 2",
      "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9),
      str_c("d", 1:13)
    ), " * irg_fraction * VPD")
  ) %>% str_flatten(collapse = "+"),
  "| x_y[year]"
) %>%
  as.formula()

fml_irrigation_heterogeneity_rtsif <- str_c(
  "RTSIF_sum ~ ",
  c(
    str_c(c(
      "Faod", "Faod ^ 2",
      "Faod:AOT40", "I(Faod ^ 2):AOT40",
      "Caod", "Caod ^ 2",
      "AOT40", "cloud", "cloud ^ 2",
      "co2",
      str_c("bin", 1:42),
      str_c("surface_", 1:9),
      str_c("d", 1:13)
    ), " * irg_fraction * VPD")
  ) %>% str_flatten(collapse = "+"),
  "| x_y[year]"
) %>%
  as.formula()