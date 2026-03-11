###########################################################
### SENSITIVITY ANALYSIS: CONTINUOUS LONELINESS SCORES ###
##########################################################

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

ESS_EULS_dataset <- read_csv("/Datasets/Data/ESS_EULS_dataset.csv")


loneliness_means <- ESS_EULS_dataset  %>% 
  dplyr::select(dataset, country, id, weight, ucla_loneli, djg_loneli, conti_dir_loneli) %>% 
  pivot_longer(-c(dataset, country, id, weight)) %>% 
  group_by(dataset, country, name) %>%
  summarise(mean = formatC(round(weighted.mean(value, w = weight, na.rm = T), 1), format = "f", digits = 1),
            sd = formatC(round((Hmisc::wtd.var(value, weights = weight, na.rm = T)/sqrt(Hmisc::wtd.var(value, weights = weight, na.rm = T))), 1), format = "f", digits = 1)) %>% 
  mutate(stat = paste0(mean, " (", sd, ")")) %>%
  mutate(stat = ifelse(stat == "NaN (NaN)", NA, stat)) %>%
  select(dataset, country, name, stat) %>% 
  pivot_wider(names_from = dataset, values_from = stat) %>% 
  mutate(value = NA) %>% 
  select(country, name, ESS11, EU27, EU4) %>% 
  arrange(name, country)

### MULTIPLE IMPUTATION ###

ESS_EULS_dataset <- ESS_EULS_dataset %>% 
  filter(dataset != "ESS7") %>% 
mutate(education = factor(education, ordered = T))

init_imp <- mice(ESS_EULS_dataset, maxit = 0) # dry run for specifying predictors and imputation model

pred_matrix <- init_imp$predictorMatrix
pred_matrix[c("id","weight", "psu", "stratum"),] <- 0 # not imputing IDs and post-stratification weights
pred_matrix[, c("id","weight", "psu", "stratum")] <- 0

imp_method <- init_imp$method
set.seed(123321)
imp_data <- mice(ESS_EULS_dataset,
                 m = 3, # number of imputed datasets
                 maxit = 3, # number of iterations
                 method = imp_method,
                 predictorMatrix = pred_matrix)

plot(imp_data)
imp_data_long_sensitivity <- complete(imp_data, "long", include = T)
write.csv(imp_data_long_sensitivity, "/Datasets/Data/imp_data_long_sensitivity.csv")




# PREVALENCE ESTIMATES 
# by data collection (ESS11, EU-LS27, EU-LS4) 
prevalence_lonely_sens <- imp_data_long_sensitivity %>%
  filter(.imp != 0) %>% 
  pivot_longer(cols = c(conti_dir_loneli, ucla_loneli, djg_loneli),
               names_to = "scale",
               values_to = "score") %>% 
  group_by(dataset, country, scale, .imp) %>% 
  nest() %>% 
  filter(dataset != "ESS7", !(dataset == "ESS11" & scale %in% c("ucla_loneli", "djg_loneli"))) %>% 
  arrange(country, dataset, .imp) %>% 
  mutate(design = map(.x = data, 
                      ~ ifelse(dataset == "ESS11", 
                               list(svydesign(ids = ~.x$psu, strata = ~.x$stratum, weights = ~.x$weight, nest = T, data = .x)),
                               list(svydesign(ids = ~.x$id, weights = ~.x$weight, data = .x)))[[1]]),
         prevalence_model = map2(.x = data,
                                 .y = design, 
                                 ~ svyglm(score ~ 1, 
                                          family = gaussian(), 
                                          data = .x,
                                          design = .y)))

pooled_prevalence_sens <- prevalence_lonely_sens %>% 
  select(prevalence_model) %>% 
  group_by(dataset, country, scale) %>% 
  nest() %>% 
  mutate(pooled = map(.x = data,
                      ~ summary(mice::pool(.x$prevalence_model), conf.int = T) %>% 
                        select(estimate, '2.5 %',  '97.5 %'))) %>% 
  select(pooled) %>% 
  unnest()
write.csv(pooled_prevalence_sens, "/Results/Tables/pooled_prevalence_sens.csv")


table_pooled_prevalence_sens <- pooled_prevalence_sens %>% 
  mutate(across(c(estimate, `2.5 %`, `97.5 %`),
                ~ formatC(round(.x, 2), format = "f", digits = 2)),
         estimate = paste0(estimate, " (", `2.5 %`, ", ", `97.5 %`, ")")) %>% 
  select(dataset, country, scale, estimate) %>% 
  pivot_wider(names_from = dataset,
              values_from = estimate)
write.csv(table_pooled_prevalence_sens, "/Results/Tables/table_pooled_prevalence_sens.csv")

windowsFonts(Times = windowsFont("Times New Roman"))
country_grey_preval <-  unique(pooled_prevalence_sens$country)[seq(1, length(unique(pooled_prevalence_sens$country)), by = 2)]

plot_prevalence_2022_conti_dir <- pooled_prevalence_sens %>% 
  filter(scale == "conti_dir_loneli") %>% 
  ungroup() %>% 
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
  labs(title = "Loneliness levels in 2022/23 by survey (sampling strategy, survey mode)",
       x = "Mean, range 1 to 5 (95% CI)",
       y = "",
       color = "") 

ggsave(filename = "plot_prevalence_2022_conti_dir.jpeg",
       path = "/Results/Graphs", 
       width = 16,
       height = 18,
       device = "jpeg",
       dpi=300)

trend_grey_country_eu4 <- c("France", "Poland")   

plot_prevalence_2022_indirect <- pooled_prevalence_sens %>% 
  filter(scale != "conti_dir_loneli", 
         country %in% c("France", "Poland", "Sweden", "Italy")) %>% 
  ungroup() %>% 
  mutate(dataset = case_when(dataset == "ESS11" ~ "ESS11, 2023 (probability sample, face-to-face interview)",
                             dataset == "EU4" ~ "EU-LS4, 2022 (probability-based panel sample, online self-completion)",
                             dataset == "EU27" ~ "EU-LS27, 2022 (nonprobability panel sample, online self-completion)"),
         scale = case_when(scale == "ucla_loneli" ~ "3-item UCLA scale", 
                           scale == "djg_loneli" ~ "6-item De Jong Gierveld scale")) %>% 
  arrange(estimate) %>% 
  ggplot(aes(x = estimate, y = country, color = dataset)) +
  geom_hline(yintercept = trend_grey_country_eu4, color = "lightgrey", size = 28, alpha = 0.5) +
  geom_point(size = 5, position = position_dodge(width = 0.75, preserve = "single")) +
  geom_errorbar(aes(xmin = `2.5 %`,
                    xmax = `97.5 %`),
                linewidth = 1,
                position = position_dodge(width = 0.75, preserve = "single"),
                width = 0.8) +
  theme_classic() +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = palette.colors(palette = "Okabe-Ito")[c(6,7)]) +  # colorblind-friendly
  facet_wrap(~ scale, 
             ncol = 1,
             scales = "free") +
  theme(legend.position = "top", legend.justification = "center",
        strip.text = element_text(size = 22, family = "Times"),
        panel.spacing = unit(5, "lines"),
        axis.text = element_text(size = 22, family = "Times"),
        axis.title.x = element_text(size = 22, family = "Times"),
        axis.text.y = element_markdown(size = 22, family = "Times", color = "black"),
        plot.title = element_text(size = 26, face = "bold", family = "Times", hjust = 0.5, margin = margin(t = 10, b = 40)),
        plot.subtitle = element_text(size = 22, family = "Times", face = "bold", hjust = 0.5, margin = margin(t = 0, b = 0)),
        legend.text = element_text(size = 22, family = "Times"), legend.margin = margin(t = 5, b = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  labs(title = "Loneliness levels in 2022/23 by survey (sampling strategy, survey mode)",
       x = "Mean, range 1 to 100 (95% CI)",
       y = "",
       color = "") 

ggsave(filename = "plot_prevalence_2022_indirect.jpeg",
       path = "/Results/Graphs", 
       width = 16,
       height = 12,
       device = "jpeg",
       dpi=300)




####### IPTW-ADJUSTMENT #########
# risk of loneliness in 2022/23 by data collection 


################################################################
### CONTINUOUS DIRECT MEASURES OF LONELINESS (I feel lonely) ###
################################################################

europe_20 <- c("Sweden", "Spain", "Slovenia", "Slovakia", "Portugal", "Poland", # countries with 2+ data collections
               "Netherlands", "Lithuania", "Italy", "Ireland", "Hungary", "Greece", 
               "Germany", "France", "Finland", "Cyprus", "Croatia", "Bulgaria", "Belgium", "Austria")

### EU-LS27 vs ESS11 ###

models_adj_EU27vsESS11_dir <- imp_data_long_sensitivity %>%
  filter(.imp != 0, 
         age >= 16 & age <= 80, 
         dataset %in% c("ESS11", "EU27"), 
         country %in% europe_20) %>% 
  mutate(dataset = relevel(as.factor(dataset), ref = "ESS11")) %>% # reference data collection, compute OR for the other datasets
  group_by(country, .imp) %>% 
  nest() %>% 
  arrange(country, .imp) %>% 
  mutate(weights = map(.x = data,
                       ~ weightit(dataset ~ age + gender + education + household_composition + marital_status + work_status + household_income + urbanicity + phys_activity + smoking + BMI,
                                  data = .x,
                                  weights = .x$weight,
                                  method = "glm",
                                  estimand = "ATE")),
         weighted_models = map2(.x = data, 
                                .y = weights,
                                ~ lm_weightit(conti_dir_loneli ~ dataset,
                                              data = .x,
                                              weightit = .y)),
         diff = map(.x = weighted_models, 
                   ~ avg_comparisons(.x, 
                                    variables = "dataset", 
                                    type = "response"))) 



pooled_adj_diff_EU27vsESS11 <- models_adj_EU27vsESS11_dir %>% 
  select(diff) %>% 
  group_by(country) %>% 
  nest() %>% 
  mutate(pooled = map(.x = data,
                      ~ summary(mice::pool(.x$diff), conf.int = T) %>%
                        filter(term != "(Intercept)") %>% 
                        select(term, estimate, '2.5 %',  '97.5 %'))) %>% 
  select(pooled) %>% 
  unnest() %>% 
  mutate(term = "EU27vsESS11")



### EU-LS4 vs ESS11 ###

models_adj_EU4vsESS11_dir <- imp_data_long_sensitivity %>%
  filter(.imp != 0, 
         age >= 16 & age <= 80,
         dataset %in% c("ESS11", "EU4"), 
         country %in% c("Sweden", "Poland",  "Italy", "France")) %>% 
  mutate(dataset = relevel(as.factor(dataset), ref = "ESS11")) %>% # reference data collection, compute OR for the other datasets
  group_by(country, .imp) %>% 
  nest() %>% 
  arrange(country, .imp) %>% 
  mutate(weights = map(.x = data,
                       ~ weightit(dataset ~ age + gender + education + household_composition + marital_status + work_status + household_income + urbanicity + phys_activity + smoking,
                                  data = .x,
                                  weights = .x$weight,
                                  method = "glm",
                                  estimand = "ATE")),
         weighted_models = map2(.x = data, 
                                .y = weights,
                                ~ lm_weightit(conti_dir_loneli ~ dataset,
                                               data = .x,
                                               weightit = .y)),
         diff = map(.x = weighted_models, 
                  ~ avg_comparisons(.x, 
                                    variables = "dataset", 
                                    type = "response")))

pooled_adj_diff_EU4vsESS11 <- models_adj_EU4vsESS11_dir %>% 
  select(diff) %>% 
  group_by(country) %>% 
  nest() %>% 
  mutate(pooled = map(.x = data,
                      ~ summary(mice::pool(.x$diff), conf.int = T) %>%
                        filter(term != "(Intercept)") %>% 
                        select(term, estimate, '2.5 %',  '97.5 %'))) %>% 
  select(pooled) %>% 
  unnest() %>% 
  mutate(term = "EU4vsESS11")




### EU-LS27 vs EU-LS4 ###

models_adj_EU27vsEU4_dir <- imp_data_long_sensitivity %>%
  filter(.imp != 0, 
         age >= 16 & age <= 80,
         dataset %in% c("EU27", "EU4"), 
         country %in% c("Sweden", "Poland",  "Italy", "France")) %>% 
  mutate(dataset = relevel(as.factor(dataset), ref = "EU4")) %>% 
  group_by(country, .imp) %>% 
  nest() %>% 
  arrange(country, .imp) %>% 
  mutate(weights = map(.x = data,
                       ~ weightit(dataset ~ age + gender + education + household_composition + marital_status + work_status + household_income + urbanicity + phys_activity + smoking,
                                  data = .x,
                                  weights = .x$weight,
                                  method = "glm",
                                  estimand = "ATE")),
         weighted_models = map2(.x = data, 
                                .y = weights,
                                ~ lm_weightit(conti_dir_loneli  ~ dataset,
                                              data = .x,
                                              weightit = .y)),
         diff = map(.x = weighted_models, 
                    ~ avg_comparisons(.x, 
                                    variables = "dataset", 
                                    type = "response")))

pooled_adj_diff_EU27vsEU4 <- models_adj_EU27vsEU4_dir %>% 
  select(diff) %>% 
  group_by(country) %>% 
  nest() %>% 
  mutate(pooled = map(.x = data,
                      ~ summary(mice::pool(.x$diff), conf.int = T) %>%
                        filter(term != "(Intercept)") %>% 
                        select(term, estimate, '2.5 %',  '97.5 %'))) %>% 
  select(pooled) %>% 
  unnest() %>% 
  mutate(term = "EU27vsEU4")



### Mini Meta-Analysis ###
# pooling effect estimates across countries into one overall estimate

all_diff_adj <- rbind(pooled_adj_diff_EU27vsESS11, pooled_adj_diff_EU4vsESS11, pooled_adj_diff_EU27vsEU4) 

meta_diff_EU4vsESS11 <- metagen(TE = estimate, 
                                lower = `2.5 %`, 
                                upper = `97.5 %`,
                                studlab = country,
                                data = pooled_adj_diff_EU4vsESS11,
                                sm = "MD")

meta_diff_EU27vsESS11 <- metagen(TE = estimate, 
                                 lower = `2.5 %`, 
                                 upper = `97.5 %`,
                                 studlab = country,
                                 data = pooled_adj_diff_EU27vsESS11,
                                 sm = "MD")

meta_diff_EU27vsEU4 <- metagen(TE = estimate, 
                               lower = `2.5 %`, 
                               upper = `97.5 %`,
                               studlab = country,
                               data = pooled_adj_diff_EU27vsEU4,
                               sm = "MD")

meta_regr_est_diff <- tibble(country = "**POOLED<br>ESTIMATE**",
                        term = c("EU4vsESS11", "EU27vsESS11", "EU27vsEU4"),
                        estimate = c(meta_diff_EU4vsESS11$TE.random, meta_diff_EU27vsESS11$TE.random, meta_diff_EU27vsEU4$TE.random), 
                        `2.5 %` = c(meta_diff_EU4vsESS11$lower.random, meta_diff_EU27vsESS11$lower.random, meta_diff_EU27vsEU4$lower.random),
                        `97.5 %` = c(meta_diff_EU4vsESS11$upper.random, meta_diff_EU27vsESS11$upper.random, meta_diff_EU27vsEU4$upper.random))

pooled_diff_conti_dir <- bind_rows(all_diff_adj, meta_regr_est_diff)
write.csv(pooled_diff_conti_dir, "/Results/Tables/pooled_diff_conti_dir.csv")

table_pooled_diff_conti_dir <- pooled_diff_conti_dir %>% 
  mutate(across(c(estimate, `2.5 %`, `97.5 %`),
                ~ formatC(round(.x, 2), format = "f", digits = 2)),
         estimate = paste0(estimate, " (", `2.5 %`, ", ", `97.5 %`, ")")) %>% 
  select(term, country, estimate) %>% 
  pivot_wider(names_from = term,
              values_from = estimate) 
write.csv(table_pooled_diff_conti_dir, "/Results/Tables/table_pooled_diff_conti_dir.csv")

windowsFonts(Times = windowsFont("Times New Roman"))
grey_country_ess11 <- c("Spain", "Slovakia", "Poland", "Lithuania", "Ireland", "Greece", "France", "Cyprus", "Bulgaria", "Austria")


plot_diff_2022_ESS11ref_conti_dir <- pooled_diff_conti_dir %>% 
  filter(str_detect(term, "vsESS11")) %>% 
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
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() +
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8)) +
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
       x = "IPTW-adjusted mean difference, range 1 to 4 (95% CI)",
       y = "",
       color = "",
       fill = "") 


plot_diff_2022_EU4ref_conti_dir <- pooled_diff_conti_dir %>% 
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
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3)) +
  theme_classic() +
  theme(legend.position = "top", legend.justification = "center", 
        axis.text = element_text(size = 22, family = "Times"),
        axis.title.x = element_text(size = 22, family = "Times"),
        axis.text.y = element_markdown(size = 22, family = "Times", color = "black"),
        plot.title = element_text(size = 26, face = "bold", family = "Times", hjust = 0.5, margin = margin(t = 10, b = 40)),
        plot.subtitle = element_text(size = 22, family = "Times", face = "bold", hjust = 0.5, margin = margin(t = 0, b = 0)),
        legend.text = element_text(size = 22, family = "Times"), legend.margin = margin(t = 5, b = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  labs(title = "Loneliness levels in 2022/23 by survey (sampling strategy, survey mode)",
       subtitle = "Reference: EU-LS4, 2022 (probability-based panel sample, online self-completion)",
       x = "IPTW-adjusted mean difference, range 1 to 4 (95% CI)",
       y = "",
       color = "",
       fill = "") 

plot_diff_2022_combined_conti_dir <- plot_diff_2022_EU4ref_conti_dir + plot_diff_2022_ESS11ref_conti_dir + 
  plot_layout(heights = c(1, 4.25))

ggsave(filename = "plot_diff_2022_combined_conti_dir.jpeg",
       path = "/Results/Graphs", 
       width = 16,
       height = 18,
       device = "jpeg",
       dpi=300)




#####################################################
### INDIRECT MULTIPLE-ITEM MEASURES OF LONELINESS ###
#####################################################


### EU-LS27 vs EU-LS4 ###

models_adj_indir <- imp_data_long_sensitivity %>%
  filter(.imp != 0, 
         age >= 16 & age <= 80,
         dataset %in% c("EU27", "EU4"), 
         country %in% c("Sweden", "Poland",  "Italy", "France")) %>% 
  pivot_longer(cols = c(conti_dir_loneli, ucla_loneli, djg_loneli),
               names_to = "scale",
               values_to = "score") %>% 
  filter(scale != "conti_dir_loneli") %>% 
  mutate(dataset = relevel(as.factor(dataset), ref = "EU4")) %>% 
  group_by(country, scale, .imp) %>% 
  nest() %>% 
  arrange(country, .imp) %>% 
  mutate(weights = map(.x = data,
                       ~ weightit(dataset ~ age + gender + education + household_composition + marital_status + work_status + household_income + urbanicity + phys_activity + smoking,
                                  data = .x,
                                  weights = .x$weight,
                                  method = "glm",
                                  estimand = "ATE")),
         weighted_models = map2(.x = data, 
                                .y = weights,
                                ~ lm_weightit(score  ~ dataset,
                                              data = .x,
                                              weightit = .y)),
         diff = map(.x = weighted_models, 
                    ~ avg_comparisons(.x, 
                                      variables = "dataset", 
                                      type = "response")))

pooled_adj_diff_indir <- models_adj_indir %>% 
  select(diff) %>% 
  group_by(country, scale) %>% 
  nest() %>% 
  mutate(pooled = map(.x = data,
                      ~ summary(mice::pool(.x$diff), conf.int = T) %>%
                        filter(term != "(Intercept)") %>% 
                        select(term, estimate, '2.5 %',  '97.5 %'))) %>% 
  select(pooled) %>% 
  unnest() %>% 
  select(-term)


### Mini Meta-Analysis ###

meta_diff_indir_ucla <- metagen(TE = estimate, 
                                lower = `2.5 %`, 
                                upper = `97.5 %`,
                                studlab = country,
                                data = pooled_adj_diff_indir[pooled_adj_diff_indir$scale == "ucla_loneli",],
                                sm = "MD")

meta_diff_indir_djg <- metagen(TE = estimate, 
                                 lower = `2.5 %`, 
                                 upper = `97.5 %`,
                                 studlab = country,
                                 data = pooled_adj_diff_indir[pooled_adj_diff_indir$scale == "djg_loneli",],
                                 sm = "MD")

meta_regr_est_diff_indir <- tibble(country = "**POOLED<br>ESTIMATE**",
                             scale = c("ucla_loneli", "djg_loneli"),
                             estimate = c(meta_diff_indir_ucla$TE.random, meta_diff_indir_djg$TE.random), 
                             `2.5 %` = c(meta_diff_indir_ucla$lower.random, meta_diff_indir_djg$lower.random),
                             `97.5 %` = c(meta_diff_indir_ucla$upper.random, meta_diff_indir_djg$upper.random))

pooled_diff_indir <- bind_rows(pooled_adj_diff_indir, meta_regr_est_diff_indir)
write.csv(pooled_diff_indir, "/Results/Tables/pooled_diff_indir.csv")

table_pooled_diff_indir <- pooled_diff_indir %>% 
  mutate(across(c(estimate, `2.5 %`, `97.5 %`),
                ~ formatC(round(.x, 2), format = "f", digits = 2)),
         estimate = paste0(estimate, " (", `2.5 %`, ", ", `97.5 %`, ")")) %>% 
  select(scale, country, estimate) %>% 
  pivot_wider(names_from = scale,
              values_from = estimate) 
write.csv(table_pooled_diff_indir, "/Results/Tables/table_pooled_diff_indir.csv")

windowsFonts(Times = windowsFont("Times New Roman"))
grey_country_ess11 <- c("Spain", "Slovakia", "Poland", "Lithuania", "Ireland", "Greece", "France", "Cyprus", "Bulgaria", "Austria")

plot_diff_indir <- pooled_diff_indir %>% 
  mutate(term = "EU-LS27, 2022 (nonprobability panel sample, online self-completion)",
         country = factor(country, levels = unique(pooled_RR_adj$country), ordered = TRUE),
         scale = case_when(scale == "ucla_loneli" ~ "3-item UCLA scale", 
                           scale == "djg_loneli" ~ "6-item De Jong Gierveld scale")) %>% 
  arrange(country) %>% 
  ggplot(aes(x = estimate, y = country, color = term, fill = term)) +
  scale_y_discrete(limits=rev) +
  geom_hline(yintercept = c("Poland", "France"), color = "lightgrey", size = 22, alpha = 0.5) +
  geom_hline(yintercept = "**POOLED<br>ESTIMATE**", color = "red", size = 22, alpha = 0.10) +
  geom_point(size = 5) +
  geom_errorbar(aes(xmin = `2.5 %`,
                    xmax = `97.5 %`),
                linewidth = 1,
                width = 0.4) +
  scale_color_manual(values = palette.colors(palette = "Okabe-Ito")[c(6,7)]) + 
  scale_fill_manual(values = palette.colors(palette = "Okabe-Ito")[c(6,7)]) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() +
  facet_wrap(~ scale, 
             ncol = 1,
             scales = "free_x") +
  theme(legend.position = "top", legend.justification = "center",
        strip.text = element_text(size = 22, family = "Times"),
        panel.spacing = unit(5, "lines"),
        axis.text = element_text(size = 22, family = "Times"),
        axis.title.x = element_text(size = 22, family = "Times"),
        axis.text.y = element_markdown(size = 22, family = "Times", color = "black"),
        plot.title = element_text(size = 26, face = "bold", family = "Times", hjust = 0.5, margin = margin(t = 10, b = 40)),
        plot.subtitle = element_text(size = 22, family = "Times", face = "bold", hjust = 0.5, margin = margin(t = 0, b = 0)),
        legend.text = element_text(size = 22, family = "Times"), legend.margin = margin(t = 5, b = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  labs(title = "Loneliness levels in 2022/23 by survey (sampling strategy, survey mode)",
       subtitle = "Reference: EU-LS4, 2022 (probability-based panel sample, online self-completion)",
       x = "IPTW-adjusted mean difference, range 1 to 100 (95% CI)",
       y = "",
       color = "",
       fill = "") 

ggsave(filename = "plot_diff_indir.jpeg",
       path = "/Results/Graphs", 
       width = 16,
       height = 12,
       device = "jpeg",
       dpi=300)






