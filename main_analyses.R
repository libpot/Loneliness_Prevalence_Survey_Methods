###########################################################################################################
### THE ROLE OF SURVEY METHODOLOGY FOR EXAMINING PREVALENCE AND TIME TRENDS IN LONELINESS ACROSS EUROPE ###
###########################################################################################################

# loading necessary packages
library(readr) # data loading
library(haven) # data loading
library(tidyverse) # data pre-processing
library(mice) # multiple imputation
library(survey) # complex survey designs
library(WeightIt) # inverse probability treatment weighting (IPTW)
library(cobalt) # covariate balance assessment
library(marginaleffects) # calculating marginal, average effects
library(sf) # visualizing geodata, maps
library(patchwork) # plot assembling
library(rnaturalearth)   # visualizing geodata, maps
library(rnaturalearthdata)  # visualizing geodata, maps
library(survey) # complex survey designs
library(meta) # meta-analysis
library(ggtext) # ggplot labels
library(modelsummary)

# loading the data
imp_data_long <- read_csv("Data/imp_data_long.csv")

# loading summary data
pooled_prevalence <- read.csv("pooled_prevalence.csv")
pooled_RR_adj <- read.csv("Tables/pooled_RR_adj.csv")
pooled_RR_adj_2023 <- read.csv("pooled_RR_adj_2023.csv")





#######################################
### THE ROLE OF SURVEY METHODOLOGY ###
######################################

# PREVALENCE ESTIMATES 
# by data collection (ESS7, ESS11, EU-LS27, EU-LS4) 
prevalence_lonely <- imp_data_long %>%
  filter(.imp != 0) %>% 
  group_by(dataset, country, .imp) %>% 
  nest() %>% 
  arrange(country, dataset, .imp) %>% 
  mutate(design = map(.x = data, 
                      ~ ifelse(dataset %in% c("ESS11", "ESS7"), 
                               list(svydesign(ids = ~.x$psu, strata = ~.x$stratum, weights = ~.x$weight, nest = T, data = .x)),
                               list(svydesign(ids = ~.x$id, weights = ~.x$weight, data = .x)))[[1]]),
         prevalence_model = map2(.x = data,
                               .y = design, 
                                ~ svyglm(I(loneliness == 1) ~ 1, 
                                         family = quasibinomial(), 
                                         data = .x,
                                         design = .y)))

pooled_prevalence <- prevalence_lonely %>% 
  select(prevalence_model) %>% 
  group_by(dataset, country) %>% 
  nest() %>% 
  mutate(pooled = map(.x = data,
                      ~ summary(mice::pool(.x$prevalence_model), conf.int = T) %>% 
                        mutate(across(c(estimate, '2.5 %',  '97.5 %'), ~ plogis(.x)*100)) %>% 
                        select(estimate, '2.5 %',  '97.5 %'))) %>% 
  select(pooled) %>% 
  unnest()
write.csv(pooled_prevalence, "/Results/Tables/pooled_prevalence.csv")


table_pooled_prevalence <- pooled_prevalence %>% 
  mutate(across(c(estimate, `2.5 %`, `97.5 %`),
                ~ formatC(round(.x, 2), format = "f", digits = 2)),
         estimate = paste0(estimate, " (", `2.5 %`, ", ", `97.5 %`, ")")) %>% 
  select(dataset, country, estimate) %>% 
  pivot_wider(names_from = dataset,
              values_from = estimate)
write.csv(table_pooled_prevalence, "/Results/Tables/table_pooled_prevalence.csv")

windowsFonts(Times = windowsFont("Times New Roman"))
country_grey_preval <-  unique(pooled_prevalence$country)[seq(1, length(unique(pooled_prevalence$country)), by = 2)]

plot_prevalence_2022 <- pooled_prevalence %>% 
  ungroup() %>% 
  filter(dataset != "ESS7") %>% 
  mutate(dataset = case_when(dataset == "ESS11" ~ "ESS11, 2023 (probability sample, face-to-face interview)",
                             dataset == "EU4" ~ "EU-LS4, 2022 (probability-based panel sample, online self-completion)",
                             dataset == "EU27" ~ "EU-LS27, 2022 (nonprobability panel sample, online self-completion)"),
         country = factor(country, ordered = TRUE)) %>% 
  arrange(estimate) %>% 
  ggplot(aes(x = estimate, y = country, color = dataset)) +
  geom_hline(yintercept = country_grey_preval, color = "lightgrey", size = 14, alpha = 0.5) +
  geom_point(size = 5, position = position_dodge(width = 0.75, preserve = "single")) +
  geom_errorbar(aes(xmin = `2.5 %`,
                    xmax = `97.5 %`),
                linewidth = 1,
                position = position_dodge(width = 0.75, preserve = "single"),
                width = 1.2) +
  theme_classic() +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = palette.colors(palette = "Okabe-Ito")[c(1,6,7)]) +  # colorblind-friendly
  theme(legend.position = "top", legend.justification = "center", 
        axis.text = element_text(size = 22, family = "Times"),
        axis.title.x = element_text(size = 22, family = "Times"),
        axis.text.y = element_markdown(size = 22, family = "Times", color = "black"),
        plot.title = element_text(size = 26, face = "bold", family = "Times", hjust = 0.5, margin = margin(t = 10, b = 40)),
        plot.subtitle = element_text(size = 22, family = "Times", face = "bold", hjust = 0.5, margin = margin(t = 0, b = 0)),
        legend.text = element_text(size = 22, family = "Times"), legend.margin = margin(t = 5, b = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  labs(title = "Loneliness prevalence in 2022/23 by survey (sampling strategy, survey mode)",
       x = "Prevalence rate (95% CI)",
       y = "",
       color = "") 

ggsave(filename = "plot_prevalence_2022.jpeg",
       path = "/Results/Graphs", 
       width = 16,
       height = 18,
       device = "jpeg",
       dpi=300)



### PREVALENCE MAP ###
pooled_prevalence_map <- select(pooled_prevalence, dataset, country, estimate) %>% 
  pivot_wider(names_from = "dataset",
              values_from = "estimate")

spatial_data <- ne_countries(scale = "medium", returnclass = "sf", continent = 'europe') %>% 
  mutate(area = st_area(geometry),
         label_point = st_point_on_surface(geometry)) %>%
  filter(area > units::set_units(2500, "km^2")) # removing small island

diff_spatial_data <- merge(spatial_data, pooled_prevalence_map, by.x = "name", by.y = "country", all.x = TRUE) 

diff_spatial_data$label_point[diff_spatial_data$name == "Norway"] <- 
  st_sfc(st_point(c(10, 61)), crs = st_crs(diff_spatial_data)) 

windowsFonts(Times = windowsFont("Times New Roman"))

ESS11_prevalence_map <- ggplot(data = diff_spatial_data) +
  geom_sf(aes(fill = ESS11)) + 
  geom_sf_text(aes(geometry = label_point,
                   label = round(ESS11, 1)),
               size = 6.5, family = "Times", color = "black") +
  scale_fill_gradient(low = palette.colors(palette = "Okabe-Ito")[4], 
                      high = palette.colors(palette = "Okabe-Ito")[7], 
                      na.value = "lightgrey",
                      limits = c(-3, 20)) + 
  labs(title = "ESS11, 2023 (probability sample, face-to-face interview)",
       fill = "") +
  theme_void() +
  coord_sf(xlim = c(-10, 34.5), ylim = c(36, 70)) +
  theme(legend.position = "none",
        plot.margin = margin(l = 25),
        plot.title = element_text(size = 22, 
                                  family = "Times",
                                  hjust = 0.5,
                                  margin = margin(b = 10, t = 5)))

ESS7_prevalence_map <- ggplot(data = diff_spatial_data) +
  geom_sf(aes(fill = ESS7)) + 
  geom_sf_text(aes(geometry = label_point,
                   label = round(ESS7, 1)),
               size = 6.5, family = "Times", color = "black") +
  scale_fill_gradient(low = palette.colors(palette = "Okabe-Ito")[4], 
                      high = palette.colors(palette = "Okabe-Ito")[7], 
                      na.value = "lightgrey",
                      limits = c(-3, 20)) + 
  labs(title = "ESS7, 2014 (probability sample, face-to-face interview)",
       fill = "") +
  theme_void() +
  coord_sf(xlim = c(-10, 34.5), ylim = c(36, 70)) +
  theme(legend.position = "none",
        plot.margin = margin(l = 25),
        plot.title = element_text(size = 22, 
                                  family = "Times",
                                  hjust = 0.5,
                                  margin = margin(b = 10, t = 5)))

EU27_prevalence_map <- ggplot(data = diff_spatial_data) +
  geom_sf(aes(fill = EU27)) + 
  geom_sf_text(aes(geometry = label_point,
                   label = round(EU27, 1)),
               size = 6.5, family = "Times", color = "black") +
  scale_fill_gradient(low = palette.colors(palette = "Okabe-Ito")[4], 
                      high = palette.colors(palette = "Okabe-Ito")[7], 
                      na.value = "lightgrey",
                      limits = c(-3, 20)) + 
  labs(title = "EU-LS27, 2022 (nonprobability panel sample, online self-completion)",
       fill = "") +
  theme_void() +
  coord_sf(xlim = c(-10, 34.5), ylim = c(36, 70)) +
  theme(legend.position = "none",
        plot.margin = margin(l = 25),
        plot.title = element_text(size = 22, 
                                  family = "Times",
                                  hjust = 0.5,
                                  margin = margin(b = 10, t = 5)))

EU4_prevalence_map <- ggplot(data = diff_spatial_data) +
  geom_sf(aes(fill = EU4)) + 
  geom_sf_text(aes(geometry = label_point,
                   label = round(EU4, 1)),
               size = 6.5, family = "Times", color = "black") +
  scale_fill_gradient(low = palette.colors(palette = "Okabe-Ito")[4], 
                      high = palette.colors(palette = "Okabe-Ito")[7], 
                      na.value = "lightgrey",
                      limits = c(-3, 20)) + 
  labs(title = "EU-LS4, 2022 (probability-based panel sample, online self-completion)",
       fill = "") +
  theme_void() +
  coord_sf(xlim = c(-10, 34.5), ylim = c(36, 70)) +
  theme(legend.position = "none",
        plot.margin = margin(l = 25),
        plot.title = element_text(size = 22, 
                                  family = "Times",
                                  hjust = 0.5,
                                  margin = margin(b = 10, t = 5)))

prevalence_map1 <- EU27_prevalence_map + ESS11_prevalence_map 
prevalence_map2 <- EU4_prevalence_map + ESS7_prevalence_map 
prevalence_map <- prevalence_map1 / prevalence_map2

ggsave(filename = "prevalence_map.jpeg",
       path = "/Results/Graphs", 
       width = 18.5,
       height = 24,
       device = "jpeg",
       dpi=400)




#################################
####### IPTW-ADJUSTMENT #########
#################################
# risk of loneliness in 2022/23 by data collection 


europe_20 <- c("Sweden", "Spain", "Slovenia", "Slovakia", "Portugal", "Poland", # countries with 2+ data collections
               "Netherlands", "Lithuania", "Italy", "Ireland", "Hungary", "Greece", 
               "Germany", "France", "Finland", "Cyprus", "Croatia", "Bulgaria", "Belgium", "Austria")

### EU-LS27 vs ESS11 ###

models_adj_EU27vsESS11 <- imp_data_long %>%
  filter(.imp != 0, 
         age >= 16 & age <= 80, 
         dataset %in% c("ESS11", "EU27"), 
         country %in% europe_20) %>% 
  mutate(dataset = relevel(as.factor(dataset), ref = "ESS11")) %>% # reference data collection, compute OR for the other datasets
  group_by(country, .imp) %>% 
  nest() %>% 
  arrange(country, .imp) %>% 
  mutate(balance_bef = map(.x = data,
                           ~ bal.tab(dataset ~ age + gender + education + household_composition + marital_status +  work_status + household_income + urbanicity + phys_activity + smoking + BMI,
                                     data = .x, 
                                     weights = .x$weight,
                                     stats = "mean.diffs",
                                     thresholds = c(m = .1))), # standardized mean differences for continuous variables and difference in proportions for binary variables
         n_nonbalanced_bef = map_int(.x = balance_bef,
                                   ~ .x$Balanced.mean.diffs["Not Balanced, >0.1",]),
         inbalance_stats_bef = map(.x = balance_bef,
                                    ~ abs(.x$Balance$Diff.Un)),
         weights = map(.x = data,
                       ~ weightit(dataset ~ age + gender + education + household_composition + marital_status + work_status + household_income + urbanicity + phys_activity + smoking + BMI,
                                  data = .x,
                                  weights = .x$weight,
                                  method = "glm",
                                  estimand = "ATE")),
         weights_trim = map(.x = weights, 
                            ~ trim(.x, at = .99)),
         prop_score_plot = map2(.x = weights,
                                .y = country,
                                ~ bal.plot(.x, 
                                           var.name = "prop.score",
                                           which = "both",
                                           type = "density",
                                           alpha = 0.8) + 
                                  theme(plot.title = element_text(hjust = 0.5),
                                        axis.text.y = element_blank(),
                                        axis.ticks.y = ggplot2::element_blank(),
                                        legend.position = "none") +
                                  labs(title = .y,
                                       y = "",
                                       x = "")),
         love_plot = map2(.x = weights,
                          .y = country,
                              ~ love.plot(.x, 
                                          binary = "std",
                                          thresholds = c(m = .1)) +
                                theme(axis.text.y = element_blank(),
                                      axis.ticks.y = ggplot2::element_blank(),
                                      legend.position = "none") +
                                labs(title = .y,
                                     x = "")),
         balance_aft = map(.x = weights,
                           ~ bal.tab(.x, stats = "mean.diffs", thresholds = c(m = .1))),
         n_nonbalanced_aft = map_int(.x = balance_aft,
                                 ~ .x$Balanced.mean.diffs["Not Balanced, >0.1",]),
         inbalance_stats_aft = map(.x = balance_aft,
                                   ~ abs(.x$Balance$Diff.Adj[-1])),
         prop_score_diff = map(.x = balance_aft, 
                              ~ as.data.frame(.x$Balance)["prop.score","Diff.Adj"]),
         data = map2(.x = data,
                     .y = weights_trim,
                    ~ .x %>% mutate(iptw = .y$weights)),
         weighted_models = map2(.x = data, 
                               .y = weights,
                               ~ glm_weightit(I(loneliness == 1) ~ dataset,
                                             data = .x,
                                             family = binomial(link = "logit"),
                                             weightit = .y)),
         RR = map(.x = weighted_models, 
                  ~ avg_comparisons(.x, 
                                    variables = "dataset", 
                                    type = "response", 
                                    comparison = "ratio"))) # risk ratios (RR)

love_plot_combined_EU27vsESS11 <- wrap_plots(filter(models_adj_EU27vsESS11, .imp == 1)$love_plot, ncol = 4) +
  plot_annotation(title = "Covariate balance before (in red) and after (in blue) IPTW-adjustment",
                  subtitle = "EU-LS27 (nonprobability, online) vs. ESS11 (probability-based, face-to-face)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                                plot.subtitle = element_text(hjust = 0.5, 
                                                             size = 13,
                                                             margin = margin(b = 10, t = 5)),
                                plot.margin = margin(5, 5, 5, 5),
                                plot.caption = element_text(hjust = 0.5, size = 13)),
                  caption = "Standardized mean difference")

ggsave(filename = "love_plot_combined_EU27vsESS11.jpeg",
       path = "/Results/Graphs", 
       width = 18,
       height = 12,
       device = "jpeg",
       dpi=300)

prop_score_plot_combined_EU27vsESS11 <- wrap_plots(filter(models_adj_EU27vsESS11, .imp == 1)$prop_score_plot, ncol = 4) +
  plot_annotation(title = "Propensity score distributions before and after IPTW-adjustment",
                  subtitle = "EU-LS27 ([in blue], agency panel, online) vs. ESS11 ([in red], probability-based, face-to-face)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                                plot.subtitle = element_text(hjust = 0.5, 
                                                             size = 13,
                                                             margin = margin(b = 10, t = 5)),
                                plot.margin = margin(5, 5, 5, 5),
                                plot.caption = element_text(hjust = 0.5, size = 13)),
                  caption = "Propensity score")

ggsave(filename = "prop_score_plot_combined_EU27vsESS11.jpeg",
       path = "/Results/Graphs", 
       width = 18,
       height = 12,
       device = "jpeg",
       dpi=300)

pooled_adj_RR_EU27vsESS11 <- models_adj_EU27vsESS11 %>% 
  select(RR) %>% 
  group_by(country) %>% 
  nest() %>% 
  mutate(pooled = map(.x = data,
                      ~ summary(mice::pool(.x$RR), conf.int = T) %>%
                        filter(term != "(Intercept)") %>% 
                        select(term, estimate, '2.5 %',  '97.5 %'))) %>% 
  select(pooled) %>% 
  unnest() %>% 
  mutate(term = "EU27vsESS11")
  


### EU-LS4 vs ESS11 ###

models_adj_EU4vsESS11 <- imp_data_long %>%
  filter(.imp != 0, 
         age >= 16 & age <= 80,
         dataset %in% c("ESS11", "EU4"), 
         country %in% c("Sweden", "Poland",  "Italy", "France")) %>% 
  mutate(dataset = relevel(as.factor(dataset), ref = "ESS11")) %>% # reference data collection, compute OR for the other datasets
  group_by(country, .imp) %>% 
  nest() %>% 
  arrange(country, .imp) %>% 
  mutate(balance_bef = map(.x = data,
                           ~ bal.tab(dataset ~ age + gender + education + household_composition + marital_status +work_status + household_income + urbanicity + phys_activity + smoking + BMI,
                                     data = .x, 
                                     weights = .x$weight,
                                     stats = "mean.diffs",
                                     thresholds = c(m = .1))), # standardized mean differences for continuous variables and difference in proportions for binary variables
         n_nonbalanced_bef = map_int(.x = balance_bef,
                                     ~ .x$Balanced.mean.diffs["Not Balanced, >0.1",]),
         inbalance_stats_bef = map(.x = balance_bef,
                                   ~ abs(.x$Balance$Diff.Un)),
         weights = map(.x = data,
                       ~ weightit(dataset ~ age + gender + education + household_composition + marital_status + work_status + household_income + urbanicity + phys_activity + smoking,
                                  data = .x,
                                  weights = .x$weight,
                                  method = "glm",
                                  estimand = "ATE")),
         weights_trim = map(.x = weights, 
                            ~ trim(.x, at = .99)),
         prop_score_plot = map2(.x = weights,
                                .y = country,
                                ~ bal.plot(.x, 
                                           var.name = "prop.score",
                                           which = "both",
                                           type = "density",
                                           alpha = 0.8) + 
                                  theme(plot.title = element_text(hjust = 0.5),
                                        axis.text.y = element_blank(),
                                        axis.ticks.y = ggplot2::element_blank(),
                                        legend.position = "none") +
                                  labs(title = .y,
                                       y = "",
                                       x = "")),
         love_plot = map2(.x = weights,
                          .y = country,
                          ~ love.plot(.x, 
                                      binary = "std",
                                      thresholds = c(m = .1)) +
                            theme(axis.text.y = element_blank(),
                                  axis.ticks.y = ggplot2::element_blank(),
                                  legend.position = "none") +
                            labs(title = .y,
                                 x = "")),
         balance_aft = map(.x = weights,
                           ~ bal.tab(.x, stats = "mean.diffs", thresholds = c(m = .1))),
         n_nonbalanced_aft = map_int(.x = balance_aft,
                                     ~ .x$Balanced.mean.diffs["Not Balanced, >0.1",]),
         inbalance_stats_aft = map(.x = balance_aft,
                                   ~ abs(.x$Balance$Diff.Adj[-1])),
         prop_score_diff = map(.x = balance_aft, 
                                   ~ as.data.frame(.x$Balance)["prop.score","Diff.Adj"]),
         data = map2(.x = data,
                     .y = weights_trim,
                     ~ .x %>% mutate(iptw = .y$weights)),
         weighted_models = map2(.x = data, 
                                .y = weights,
                                ~ glm_weightit(I(loneliness == 1) ~ dataset,
                                               data = .x,
                                               family = binomial(link = "logit"),
                                               weightit = .y)),
         RR = map(.x = weighted_models, 
                  ~ avg_comparisons(.x, 
                                    variables = "dataset", 
                                    type = "response", 
                                    comparison = "ratio")))


love_plot_combined_EU4vsESS11 <- wrap_plots(filter(models_adj_EU4vsESS11, .imp == 1)$love_plot, ncol = 2) +
  plot_annotation(title = "Covariate balance before (in red) and after (in blue) IPTW-adjustment",
                  subtitle = "EU-LS4 (probability-based, online) vs. ESS11 (probability-based, face-to-face)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                                plot.subtitle = element_text(hjust = 0.5, 
                                                             size = 13,
                                                             margin = margin(b = 10, t = 5)),
                                plot.margin = margin(5, 5, 5, 5),
                                plot.caption = element_text(hjust = 0.5, size = 13)),
                  caption = "Standardized mean difference")

ggsave(filename = "love_plot_combined_EU4vsESS11.jpeg",
       path = "/Results/Graphs", 
       width = 9,
       height = 9,
       device = "jpeg",
       dpi=200)

prop_score_plot_combined_EU4vsESS11 <- wrap_plots(filter(models_adj_EU4vsESS11, .imp == 1)$prop_score_plot, ncol = 2) +
  plot_annotation(title = "Propensity score distributions before and after IPTW-adjustment",
                  subtitle = "EU-LS4 ([in blue], probability-based, online) vs. ESS11 ([in red], probability-based, face-to-face)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                                plot.subtitle = element_text(hjust = 0.5, 
                                                             size = 13,
                                                             margin = margin(b = 10, t = 5)),
                                plot.margin = margin(5, 5, 5, 5),
                                plot.caption = element_text(hjust = 0.5, size = 13)),
                  caption = "Propensity score")

ggsave(filename = "prop_score_plot_combined_EU4vsESS11.jpeg",
       path = "/Results/Graphs", 
       width = 9,
       height = 9,
       device = "jpeg",
       dpi=200)

pooled_adj_RR_EU4vsESS11 <- models_adj_EU4vsESS11 %>% 
  select(RR) %>% 
  group_by(country) %>% 
  nest() %>% 
  mutate(pooled = map(.x = data,
                      ~ summary(mice::pool(.x$RR), conf.int = T) %>%
                        filter(term != "(Intercept)") %>% 
                        select(term, estimate, '2.5 %',  '97.5 %'))) %>% 
  select(pooled) %>% 
  unnest() %>% 
  mutate(term = "EU4vsESS11")




### EU-LS27 vs EU-LS4 ###

models_adj_EU27vsEU4 <- imp_data_long %>%
  filter(.imp != 0, 
         age >= 16 & age <= 80,
         dataset %in% c("EU27", "EU4"), 
         country %in% c("Sweden", "Poland",  "Italy", "France")) %>% 
  mutate(dataset = relevel(as.factor(dataset), ref = "EU4")) %>% 
  group_by(country, .imp) %>% 
  nest() %>% 
  arrange(country, .imp) %>% 
  mutate(balance_bef = map(.x = data,
                           ~ bal.tab(dataset ~ age + gender + education + household_composition + marital_status +work_status + household_income + urbanicity + phys_activity + smoking + BMI,
                                     data = .x, 
                                     weights = .x$weight,
                                     stats = "mean.diffs",
                                     thresholds = c(m = .1))), # standardized mean differences for continuous variables and difference in proportions for binary variables
         n_nonbalanced_bef = map_int(.x = balance_bef,
                                     ~ .x$Balanced.mean.diffs["Not Balanced, >0.1",]),
         inbalance_stats_bef = map(.x = balance_bef,
                                   ~ abs(.x$Balance$Diff.Un)),
         weights = map(.x = data,
                       ~ weightit(dataset ~ age + gender + education + household_composition + marital_status + work_status + household_income + urbanicity + phys_activity + smoking,
                                  data = .x,
                                  weights = .x$weight,
                                  method = "glm",
                                  estimand = "ATE")),
         weights_trim = map(.x = weights, 
                            ~ trim(.x, at = .99)),
         prop_score_plot = map2(.x = weights,
                                .y = country,
                                ~ bal.plot(.x, 
                                           var.name = "prop.score",
                                           which = "both",
                                           type = "density",
                                           alpha = 0.8) + 
                                  theme(plot.title = element_text(hjust = 0.5),
                                        axis.text.y = element_blank(),
                                        axis.ticks.y = ggplot2::element_blank(),
                                        legend.position = "none") +
                                  labs(title = .y,
                                       y = "",
                                       x = "")),
         love_plot = map2(.x = weights,
                          .y = country,
                          ~ love.plot(.x, 
                                      binary = "std",
                                      thresholds = c(m = .1)) +
                            theme(axis.text.y = element_blank(),
                                  axis.ticks.y = ggplot2::element_blank(),
                                  legend.position = "none") +
                            labs(title = .y,
                                 x = "")),
         balance_aft = map(.x = weights,
                           ~ bal.tab(.x, stats = "mean.diffs", thresholds = c(m = .1))),
         n_nonbalanced_aft = map_int(.x = balance_aft,
                                     ~ .x$Balanced.mean.diffs["Not Balanced, >0.1",]),
         inbalance_stats_aft = map(.x = balance_aft,
                                   ~ abs(.x$Balance$Diff.Adj[-1])),
         prop_score_diff = map(.x = balance_aft, 
                                   ~ as.data.frame(.x$Balance)["prop.score","Diff.Adj"]),
         data = map2(.x = data,
                     .y = weights_trim,
                     ~ .x %>% mutate(iptw = .y$weights)),
         weighted_models = map2(.x = data, 
                                .y = weights,
                                ~ glm_weightit(I(loneliness == 1) ~ dataset,
                                               data = .x,
                                               family = binomial(link = "logit"),
                                               weightit = .y)),
         RR = map(.x = weighted_models, 
                  ~ avg_comparisons(.x, 
                                    variables = "dataset", 
                                    type = "response", 
                                    comparison = "ratio")))


love_plot_combined_EU27vsEU4 <- wrap_plots(filter(models_adj_EU27vsEU4, .imp == 1)$love_plot, ncol = 2) +
  plot_annotation(title = "Covariate balance before (in red) and after (in blue) IPTW-adjustment",
                  subtitle = "EU-LS27 (nonprobability, online) vs. EU-LS4 (probability-based, online)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                                plot.subtitle = element_text(hjust = 0.5, 
                                                             size = 13,
                                                             margin = margin(b = 10, t = 5)),
                                plot.margin = margin(5, 5, 5, 5),
                                plot.caption = element_text(hjust = 0.5, size = 13)),
                  caption = "Standardized mean difference")

ggsave(filename = "love_plot_combined_EU27vsEU4.jpeg",
       path = "/Results/Graphs", 
       width = 9,
       height = 9,
       device = "jpeg",
       dpi=200)

prop_score_plot_combined_EU27vsEU4 <- wrap_plots(filter(models_adj_EU27vsEU4, .imp == 1)$prop_score_plot, ncol = 2) +
  plot_annotation(title = "Propensity score distributions before and after IPTW-adjustment",
                  subtitle = "EU-LS27 ([in blue], nonprobability, online) vs. EU-LS4 ([in red], probability-based, online)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                                plot.subtitle = element_text(hjust = 0.5, 
                                                             size = 13,
                                                             margin = margin(b = 10, t = 5)),
                                plot.margin = margin(5, 5, 5, 5),
                                plot.caption = element_text(hjust = 0.5, size = 13)),
                  caption = "Propensity score")

ggsave(filename = "prop_score_plot_combined_EU27vsEU4.jpeg",
       path = "/Results/Graphs", 
       width = 9,
       height = 9,
       device = "jpeg",
       dpi=200)


pooled_adj_RR_EU27vsEU4 <- models_adj_EU27vsEU4 %>% 
  select(RR) %>% 
  group_by(country) %>% 
  nest() %>% 
  mutate(pooled = map(.x = data,
                      ~ summary(mice::pool(.x$RR), conf.int = T) %>%
                        filter(term != "(Intercept)") %>% 
                        select(term, estimate, '2.5 %',  '97.5 %'))) %>% 
  select(pooled) %>% 
  unnest() %>% 
  mutate(term = "EU27vsEU4")



### Mini Meta-Analysis ###
# pooling effect estimates across countries into one overall estimate

all_RR_adj <- rbind(pooled_adj_RR_EU27vsESS11, pooled_adj_RR_EU4vsESS11, pooled_adj_RR_EU27vsEU4) 

meta_regr_EU4vsESS11 <- metagen(TE = log(estimate), 
                                lower = log(`2.5 %`), 
                                upper = log(`97.5 %`),
                                studlab = country,
                                data = pooled_adj_RR_EU4vsESS11,
                                sm = "RR")

meta_regr_EU27vsESS11 <- metagen(TE = log(estimate), 
                                lower = log(`2.5 %`), 
                                upper = log(`97.5 %`),
                                studlab = country,
                                data = pooled_adj_RR_EU27vsESS11,
                                sm = "RR")

meta_regr_EU27vsEU4 <- metagen(TE = log(estimate), 
                               lower = log(`2.5 %`), 
                               upper = log(`97.5 %`),
                               studlab = country,
                               data = pooled_adj_RR_EU27vsEU4,
                               sm = "RR")

meta_regr_est <- tibble(country = "**POOLED<br>ESTIMATE**",
                        term = c("EU4vsESS11", "EU27vsESS11", "EU27vsEU4"),
                        estimate = c(exp(meta_regr_EU4vsESS11$TE.random), exp(meta_regr_EU27vsESS11$TE.random), exp(meta_regr_EU27vsEU4$TE.random)), 
                        `2.5 %` = c(exp(meta_regr_EU4vsESS11$lower.random), exp(meta_regr_EU27vsESS11$lower.random), exp(meta_regr_EU27vsEU4$lower.random)),
                        `97.5 %` = c(exp(meta_regr_EU4vsESS11$upper.random), exp(meta_regr_EU27vsESS11$upper.random), exp(meta_regr_EU27vsEU4$upper.random)))

pooled_RR_adj <- bind_rows(all_RR_adj, meta_regr_est)
write.csv(pooled_RR_adj, "/Results/Tables/pooled_RR_adj.csv")

table_pooled_RR_adj <- pooled_RR_adj %>% 
  mutate(across(c(estimate, `2.5 %`, `97.5 %`),
                ~ formatC(round(.x, 2), format = "f", digits = 2)),
         estimate = paste0(estimate, " (", `2.5 %`, ", ", `97.5 %`, ")")) %>% 
  select(term, country, estimate) %>% 
  pivot_wider(names_from = term,
              values_from = estimate) 
write.csv(table_pooled_RR_adj, "/Results/Tables/table_pooled_RR_adj.csv")

windowsFonts(Times = windowsFont("Times New Roman"))
grey_country_ess11 <- c("Spain", "Slovakia", "Poland", "Lithuania", "Ireland", "Greece", "France", "Cyprus", "Bulgaria", "Austria")


plot_RR_2022_ESS11ref <- pooled_RR_adj %>% 
  filter(str_detect(term, "vsESS11")) %>% 
  mutate(`2.5 %` = `X2.5..`,
         `97.5 %` = `X97.5..`) %>% 
  mutate(term = case_when(term == "EU27vsESS11" ~ "EU-LS27, 2022 (nonprobability panel sample, online self-completion)",
                          term == "EU4vsESS11" ~ "EU-LS4, 2022 (probability-based panel sample, online self-completion)"),
         country = factor(country, levels = unique(pooled_RR_adj$country), ordered = TRUE)) %>% 
  arrange(country) %>% 
  ggplot(aes(x = estimate, y = country, color = term, fill = term)) +
  scale_y_discrete(limits=rev) +
  geom_hline(yintercept = grey_country_ess11, color = "lightgrey", size = 16, alpha = 0.5) +
  geom_hline(yintercept = "**POOLED<br>ESTIMATE**", color = "red", size = 16, alpha = 0.10) +
  geom_point(position = position_dodge(width = 0.55, preserve = "single"),
             size = 5) +
  geom_errorbar(aes(xmin = `2.5 %`,
                    xmax = `97.5 %`),
                position = position_dodge(width = 0.55, preserve = "single"),
                linewidth = 1,
                width = 0.8) +
  scale_color_manual(values = palette.colors(palette = "Okabe-Ito")[c(6,7)]) + 
  scale_fill_manual(values = palette.colors(palette = "Okabe-Ito")[c(6,7)]) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  theme_classic() +
  scale_x_continuous(breaks = c(0.5,0,1,2,3)) +
  theme(legend.position = "top", legend.justification = "center",
        plot.margin = margin(t = 40, r = 0, b = 0, l = 0),
        axis.text = element_text(size = 22, family = "Times"),
        axis.title.x = element_text(size = 22, family = "Times"),
        axis.text.y = element_markdown(size = 22, family = "Times", color = "black"),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size = 22, face = "bold", family = "Times", hjust = 0.5, margin = margin(t = 2, b = 10)),
        legend.text = element_text(size = 22, family = "Times"),
        legend.box.margin = margin(t = -5, b = -10),
        legend.margin = margin(b = 20)) +
  guides(color = guide_legend(ncol = 1)) +
  labs(subtitle = "Reference: ESS11, 2023 (probability sample, face-to-face interview)",
       x = "IPTW-adjusted risk ratio (95% CI)",
       y = "",
       color = "",
       fill = "") 


plot_RR_2022_EU4ref <- pooled_RR_adj %>% 
  mutate(`2.5 %` = `X2.5..`,
         `97.5 %` = `X97.5..`) %>% 
  filter(str_detect(term, "vsEU4")) %>% 
  mutate(term = case_when(term == "EU27vsEU4" ~ "EU-LS27, 2022 (nonprobability panel sample, online self-completion)"),
         country = factor(country, levels = unique(pooled_RR_adj$country), ordered = TRUE)) %>% 
  arrange(country) %>% 
  ggplot(aes(x = estimate, y = country, color = term, fill = term)) +
  scale_y_discrete(limits=rev) +
  geom_hline(yintercept = c("Poland", "France"), color = "lightgrey", size = 16, alpha = 0.5) +
  geom_hline(yintercept = "**POOLED<br>ESTIMATE**", color = "red", size = 16, alpha = 0.10) +
  geom_point(size = 5) +
  geom_errorbar(aes(xmin = `2.5 %`,
                    xmax = `97.5 %`),
                linewidth = 1,
                width = 0.4) +
  scale_color_manual(values = palette.colors(palette = "Okabe-Ito")[c(6,7)]) + 
  scale_fill_manual(values = palette.colors(palette = "Okabe-Ito")[c(6,7)]) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  theme_classic() +
  theme(legend.position = "top", legend.justification = "center", 
        axis.text = element_text(size = 22, family = "Times"),
        axis.title.x = element_text(size = 22, family = "Times"),
        axis.text.y = element_markdown(size = 22, family = "Times", color = "black"),
        plot.subtitle = element_text(size = 22, family = "Times", face = "bold", hjust = 0.5, margin = margin(t = 2, b = 10)),
        legend.text = element_text(size = 22, family = "Times"), legend.margin = margin(t = 5, b = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  labs(subtitle = "Reference: EU-LS4, 2022 (probability-based panel sample, online self-completion)",
       x = "IPTW-adjusted risk ratio (95% CI)",
       y = "",
       color = "",
       fill = "") 

plot_RR_2022_combined <- plot_RR_2022_EU4ref + plot_RR_2022_ESS11ref + 
  plot_layout(heights = c(1, 4.25))

ggsave(filename = "plot_RR_2022_combined.jpeg",
       path = "/Results/Graphs", 
       width = 16,
       height = 18,
       device = "jpeg",
       dpi=300)


plot_2022_combined <- plot_prevalence_2022 + plot_RR_2022_combined 

ggsave(filename = "plot_2022_combined.jpeg",
       path = "/Results/Graphs", 
       width = 30,
       height = 18,
       device = "jpeg",
       dpi=300)



  


########################
### TEMPORAL TRENDS ####
#######################

time_diff_prevalence <- pooled_prevalence %>% 
  select(dataset, country, estimate) %>% 
  pivot_wider(names_from = dataset,
              values_from = estimate) %>% 
  mutate(ess11_vs_ess7 = ESS11 - ESS7,
         eu27_vs_ess7 = EU27 - ESS7,
         eu4_vs_ess7 = EU4 - ESS7)

spatial_data <- ne_countries(scale = "medium", returnclass = "sf", continent = 'europe') %>% 
  mutate(area = st_area(geometry),
         label_point = st_point_on_surface(geometry)) %>%
  filter(area > units::set_units(2500, "km^2")) # removing small island

diff_spatial_data <- merge(spatial_data, time_diff_prevalence, by.x = "name", by.y = "country", all.x = TRUE) 

diff_spatial_data$label_point[diff_spatial_data$name == "Norway"] <- 
  st_sfc(st_point(c(10, 61)), crs = st_crs(diff_spatial_data))


windowsFonts(Times = windowsFont("Times New Roman"))

diff_plot_ess11 <- ggplot(data = diff_spatial_data) +
  geom_sf(aes(fill = ess11_vs_ess7)) +  
  geom_sf_text(aes(geometry = label_point,
                   label = round(ess11_vs_ess7, 1)),
               size = 5, family = "Times", color = "black") +
  scale_fill_gradient(low = palette.colors(palette = "Okabe-Ito")[4], 
                      high = palette.colors(palette = "Okabe-Ito")[7], 
                      na.value = "lightgrey",
                      limits = c(-3, 15), breaks = c(0, 5, 10)) + 
  labs(title = "ESS11, 2023 (probability sample, face-to-face interview)",
       fill = "") +
  theme_void() +
  coord_sf(xlim = c(-10, 34.5), ylim = c(36, 70)) +
  theme(legend.position = "bottom",
        plot.margin = margin(l = 25),
        plot.title = element_text(size = 19, 
                                  family = "Times",
                                  hjust = 0.5,
                                  margin = margin(b = 10, t = 5)),
        legend.text = element_text(size = 19, family = "Times"))

diff_plot_eu27 <- ggplot(data = diff_spatial_data) +
  geom_sf(aes(fill = eu27_vs_ess7)) +  
  geom_sf_text(aes(geometry = label_point,
                   label = round(eu27_vs_ess7, 1)),
               size = 5, family = "Times", color = "black") +
  scale_fill_gradient(low = palette.colors(palette = "Okabe-Ito")[4], 
                      high = palette.colors(palette = "Okabe-Ito")[7], 
                      na.value = "lightgrey",
                      limits = c(-3, 15), breaks = c(0, 5, 10)) +  
  labs(title = "EU-LS27, 2022 (nonprobability panel sample, online self-completion)",
       fill = "") +
  theme_void() +
  coord_sf(xlim = c(-10, 34.5), ylim = c(36, 70)) +
  theme(legend.position = "bottom",
        plot.margin = margin(r = 25),
        plot.title = element_text(size = 19, 
                                  family = "Times",
                                  hjust = 0.5,
                                  margin = margin(b = 10, t = 5)),
        legend.text = element_text(size = 19, family = "Times"))

windowsFonts(Times = windowsFont("Times New Roman"))


plot_prevalence_trend <- diff_plot_eu27 + diff_plot_ess11 +  
  plot_annotation(title = "Change in loneliness prevalence between 2014 and 2022/23",
                  subtitle = " Reference: ESS7, 2014 (probability sample, face-to-face interview)",
                  theme = theme(plot.title = element_text(size = 22, 
                                             family = "Times",
                                             face = "bold", 
                                             hjust = 0.5,
                                             margin = margin(b = 35, t = 10)),
                                plot.subtitle = element_text(size = 19, 
                                                             family = "Times",
                                                             face = "bold", 
                                                             hjust = 0.5),
                                plot.caption = element_text(hjust = 0.5, size = 19, family = "Times")),
                  labs(caption = "Change in percentage points (%)"))


ggsave(filename = "plot_prevalence_trend.jpeg",
       path = "/Results/Graphs", 
       width = 18,
       height = 12,
       device = "jpeg",
       dpi=300)



ess11_countries <- filter(time_diff_prevalence, !is.na(ess11_vs_ess7))$country
eu4_countries <- filter(time_diff_prevalence, !is.na(eu4_vs_ess7))$country
eu27_countries <- filter(time_diff_prevalence, !is.na(eu27_vs_ess7))$country

# risk of loneliness in 2022/23 (EU-LS27, EU-LS4) vs 2014 (ESS7)
RR_2023_vs_2014 <- imp_data_long %>%
  filter(.imp != 0, 
         country %in% ess11_countries | country %in% eu27_countries) %>% # countries with both 2014 and 2022/23
  #mutate(dataset = factor(dataset, levels = c("ESS7", "ESS11", "EU4", "EU27"))) %>% # reference data collection, compute OR for the other datasets
  group_by(country, dataset, .imp) %>% 
  nest() %>% 
  mutate(data = map(.x = data,
                    .y = dataset,
                    ~ mutate(.x, dataset = .y))) %>% 
  pivot_wider(names_from = dataset,
              values_from = data) %>% 
  mutate(EU27 = map2(.x = EU27,
                     .y = ESS7,
                     ~ bind_rows(.x, .y)),
         EU4 = map2(.x = EU4,
                     .y = ESS7,
                     ~ bind_rows(.x, .y)),
         ESS11 = map2(.x = ESS11,
                     .y = ESS7,
                     ~ bind_rows(.x, .y))) %>%
  select(-ESS7) %>% 
  pivot_longer(cols = ESS11:EU4,
               names_to = "dataset",
               values_to = "data") %>% 
  filter(dataset == "ESS11" & country %in% ess11_countries | 
         dataset == "EU27" & country %in% eu27_countries | 
         dataset == "EU4" & country %in% eu4_countries) %>% 
  mutate(data = map(.x = data,
                    ~ mutate(.x, exposure = ifelse(dataset == "ESS7", 0, 1))),
         design = map(.x = data, 
                      ~ ifelse(.x$dataset %in% c("ESS11", "ESS7"), 
                               list(svydesign(ids = ~.x$psu, strata = ~.x$stratum, weights = ~.x$weight, nest = T, data = .x)),
                               list(svydesign(ids = ~.x$id, weights = ~.x$weight, data = .x)))[[1]]),
         outcome_model = map2(.x = data, 
                                    .y = design, 
                                    ~ svyglm(I(loneliness == 1) ~  exposure, 
                                             family = quasibinomial, 
                                             data = .x,
                                             design = .y)),
         RR = map2(.x = outcome_model, 
                   .y = data,
                   ~ avg_comparisons(.x, 
                                     variables = "exposure", 
                                     type = "response", 
                                     comparison = "ratio",
                                     wts = .y$weight)),
         RD = map2(.x = outcome_model, 
                   .y = data,
                   ~ avg_comparisons(.x, 
                                     variables = "exposure", 
                                     type = "response", 
                                     comparison = "difference",
                                     wts = .y$weight)))

pooled_RR_2023_vs_2014 <- RR_2023_vs_2014 %>% 
  select(country, .imp, dataset, RR) %>% 
  group_by(country, dataset) %>% 
  nest() %>% 
  mutate(RR = map(.x = data,
                  ~ summary(mice::pool(.x$RR), conf.int = T) %>%
                  filter(term != "(Intercept)") %>% 
                  select(estimate, '2.5 %',  '97.5 %'))) %>% 
  select(RR) %>% 
  unnest() %>% 
  arrange(country)

pooled_RD_2023_vs_2014 <- RR_2023_vs_2014 %>% 
  select(country, .imp, dataset, RD) %>% 
  group_by(country, dataset) %>% 
  nest() %>% 
  mutate(RD = map(.x = data,
                  ~ summary(mice::pool(.x$RD), conf.int = T) %>%
                    filter(term != "(Intercept)") %>% 
                    select(estimate, '2.5 %',  '97.5 %') %>% 
                    mutate(across(c(estimate, '2.5 %',  '97.5 %'),
                                  ~ .x*100)))) %>% 
  select(RD) %>% 
  unnest() %>% 
  arrange(country)

meta_RR_2023 <- pooled_RR_2023_vs_2014 %>% 
  group_by(dataset) %>% 
  nest() %>% 
  mutate(meta_regr = map(.x = data,
                         ~ metagen(TE = log(.x$estimate), 
                                   lower = log(.x$`2.5 %`), 
                                   upper = log(.x$`97.5 %`),
                                   studlab = .x$country,
                                   data = .x,
                                   sm = "RR")),
        pooled_RR = map(.x = meta_regr,
                        ~ tibble(country = "**POOLED<br>ESTIMATE**",
                                 estimate = exp(.x$TE.random),  
                                 `2.5 %` = exp(.x$lower.random), 
                                 `97.5 %` = exp(.x$upper.random)))) %>% 
  select(pooled_RR) %>% 
  unnest()

meta_RD_2023 <- pooled_RD_2023_vs_2014 %>% 
  group_by(dataset) %>% 
  nest() %>% 
  mutate(meta_regr = map(.x = data,
                          ~ metagen(TE = .x$estimate, 
                                    lower = .x$`2.5 %`, 
                                    upper = .x$`97.5 %`,
                                    studlab = .x$country,
                                    data = .x,
                                    sm = "RD")),
         pooled_RD = map(.x = meta_regr,
                          ~ tibble(country = "**POOLED<br>ESTIMATE**",
                                   estimate = .x$TE.random,  
                                   `2.5 %` = .x$lower.random, 
                                   `97.5 %` = .x$upper.random))) %>% 
  select(pooled_RD) %>% 
  unnest()

pooled_RR_adj_2023 <- bind_rows(pooled_RR_2023_vs_2014, meta_RR_2023)
pooled_RD_adj_2023 <- bind_rows(pooled_RD_2023_vs_2014, meta_RD_2023)
write.csv(pooled_RR_adj_2023, "/Results/Tables/pooled_RR_adj_2023.csv")
write.csv(pooled_RD_adj_2023, "/Results/Tables/pooled_RD_adj_2023.csv")

table_pooled_RR_adj_2023 <- pooled_RR_adj_2023 %>% 
  mutate(across(c(estimate, `2.5 %`, `97.5 %`),
                ~ formatC(round(.x, 2), format = "f", digits = 2)),
         estimate = paste0(estimate, " (", `2.5 %`, ", ", `97.5 %`, ")")) %>% 
  select(dataset, country, estimate) %>% 
  pivot_wider(names_from = dataset,
              values_from = estimate) 
write.csv(table_pooled_RR_adj_2023, "/Results/Tables/table_pooled_RR_adj_2023.csv")

table_pooled_RD_adj_2023 <- pooled_RD_adj_2023 %>% 
  mutate(across(c(estimate, `2.5 %`, `97.5 %`),
                ~ formatC(round(.x, 2), format = "f", digits = 2)),
         estimate = paste0(estimate, " (", `2.5 %`, ", ", `97.5 %`, ")")) %>% 
  select(dataset, country, estimate) %>% 
  pivot_wider(names_from = dataset,
              values_from = estimate) 
write.csv(table_pooled_RD_adj_2023, "/Results/Tables/table_pooled_RD_adj_2023.csv")


trend_grey_country <- c("Austria", "Czechia", "Estonia", "France", "Hungary", "Israel", "Netherlands", "Poland", "Slovenia", "Sweden", "United Kingdom")   

plot_RR_2023 <- pooled_RR_adj_2023 %>% 
  mutate(`2.5 %` = `X2.5..`,
         `97.5 %` = `X97.5..`) %>% 
  arrange(country) %>% 
  mutate(dataset = case_when(dataset == "ESS11" ~ "ESS11, 2023 (probability sample, face-to-face interview)",
                             dataset == "EU27" ~ "EU-LS27, 2022 (nonprobability panel sample, online self-completion)",
                             dataset == "EU4" ~ "EU-LS4, 2022 (probability-based panel sample, online self-completion)"),
         country = factor(country, levels = unique(pooled_RR_adj_2023$country), ordered = TRUE)) %>% 
  ggplot(aes(x = estimate, y = country, color = dataset)) +
  scale_y_discrete(limits=rev) +
  geom_hline(yintercept = trend_grey_country, color = "lightgrey", size = 19, alpha = 0.5) +
  geom_hline(yintercept = "**POOLED<br>ESTIMATE**", color = "red", size = 19, alpha = 0.10) +
  geom_point(size = 5, position = position_dodge(width = 0.75, preserve = "single")) +
  geom_errorbar(aes(xmin = `2.5 %`,
                    xmax = `97.5 %`),
                linewidth = 1,
                position = position_dodge(width = 0.75, preserve = "single"),
                width = 0.8) +
  scale_color_manual(values = palette.colors(palette = "Okabe-Ito")[c(1,6,7)]) +  # colorblind-friendly
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  theme_classic() +
  theme(legend.position = "top", legend.justification = "center", 
        axis.text = element_text(size = 22, family = "Times"),
        axis.title.x = element_text(size = 22, family = "Times"),
        axis.text.y = element_markdown(size = 22, family = "Times", color = "black"),
        #plot.title = element_text(size = 26, face = "bold", family = "Times", hjust = 0.5, margin = margin(t = 10, b = 40)),
        plot.subtitle = element_text(size = 22, family = "Times", face = "bold", hjust = 0.5, margin = margin(t = 0, b = 0)),
        legend.text = element_text(size = 22, family = "Times"), legend.margin = margin(t = 5, b = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  labs(#title = "Loneliness risk in 2022/23 vs. 2014 by survey (sampling strategy, survey mode)",
       subtitle = "Reference: ESS7, 2014 (probability sample, face-to-face interview)",
       x = "Risk ratio (95% CI)",
       y = "",
       color = "") 

ggsave(filename = "plot_RR_2023.jpeg",
       path = "/Results/Graphs", 
       width = 16,
       height = 16,
       device = "jpeg",
       dpi=300)


plot_RD_2023 <- pooled_RD_adj_2023 %>% 
  arrange(country) %>% 
  mutate(dataset = case_when(dataset == "ESS11" ~ "ESS11, 2023 (probability sample, face-to-face interview)",
                             dataset == "EU27" ~ "EU-LS27, 2022 (nonprobability panel sample, online self-completion)",
                             dataset == "EU4" ~ "EU-LS4, 2022 (probability-based panel sample, online self-completion)"),
         country = factor(country, levels = unique(pooled_RR_adj_2023$country), ordered = TRUE)) %>% 
  ggplot(aes(x = estimate, y = country, color = dataset)) +
  scale_y_discrete(limits=rev) +
  geom_hline(yintercept = trend_grey_country, color = "lightgrey", size = 19, alpha = 0.5) +
  geom_hline(yintercept = "**POOLED<br>ESTIMATE**", color = "red", size = 19, alpha = 0.10) +
  geom_point(size = 5, position = position_dodge(width = 0.75, preserve = "single")) +
  geom_errorbar(aes(xmin = `2.5 %`,
                    xmax = `97.5 %`),
                linewidth = 1,
                position = position_dodge(width = 0.75, preserve = "single"),
                width = 0.8) +
  scale_color_manual(values = palette.colors(palette = "Okabe-Ito")[c(1,6,7)]) +  # colorblind-friendly
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() +
  theme(legend.position = "top", legend.justification = "center", 
        axis.text = element_text(size = 22, family = "Times"),
        axis.title.x = element_text(size = 22, family = "Times"),
        axis.text.y = element_markdown(size = 22, family = "Times", color = "black"),
        plot.title = element_text(size = 26, face = "bold", family = "Times", hjust = 0.5, margin = margin(t = 10, b = 40)),
        plot.subtitle = element_text(size = 22, family = "Times", face = "bold", hjust = 0.5, margin = margin(t = 0, b = 0)),
        legend.text = element_text(size = 22, family = "Times"), legend.margin = margin(t = 5, b = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  labs(title = "Loneliness risk in 2022/23 vs. 2014 by survey (sampling strategy, survey mode)",
       subtitle = "Reference: ESS7, 2014 (probability sample, face-to-face interview)",
       x = "Risk difference (95% CI)",
       y = "",
       color = "") 

ggsave(filename = "plot_RD_2023.jpeg",
       path = "/Results/Graphs", 
       width = 16,
       height = 16,
       device = "jpeg",
       dpi=300)

