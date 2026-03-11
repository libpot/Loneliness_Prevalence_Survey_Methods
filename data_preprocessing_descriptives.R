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
ESS11 <- read_csv("C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Datasets/Data/ESS11/ESS11.csv")
EU27 <- read_dta("C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Datasets/Data/EU_Lonely/eu_loneliness_survey_eu27.dta")
EU4 <- read_dta("C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Datasets/Data/EU_Lonely_E4/eu_loneliness_survey_eu4.dta")
ESS7 <- read_sav("C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Datasets/Data/ESS7/ESS7e02_3.sav")
ESS7_w <- read_sav("C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Datasets/Data/ESS7/ESS7SDDFe1_2.sav")
ESS7 <- as_tibble(merge(ESS7, ESS7_w, by = c("cntry", "idno")))
imp_data_long <- read_csv("C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Datasets/Data/imp_data_long.csv")

# loading summary data
pooled_prevalence <- read.csv("C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Prevalence/Results/Tables/pooled_prevalence.csv")
pooled_RR_adj <- read.csv("C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Prevalence/Results/Tables/pooled_RR_adj.csv")
pooled_RR_adj_2023 <- read.csv("C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Prevalence/Results/Tables/pooled_RR_adj_2023.csv")



###########################
### DATA PRE-PROCESSING ###
###########################
#   data harmonization   #


# European Social Survey 11 (ESS11, 2023, probability-based household sample, in-person interviewed)
ESS11_data <- ESS11 %>%                               # How much of the time during the past week...you felt lonely?
  mutate(loneliness = case_when(fltlnl %in% c(3,4) ~ 1, # all or most of time
                                fltlnl %in% c(1,2) ~ 0, # none or some of time
                                fltlnl %in% c(7,8,9) ~ NA), # refusal, don't know or no answer
         conti_dir_loneli = case_when(fltlnl %in% c(1,2,3,4) ~ fltlnl, # ~ 1 + (fltlnl - 1) * 4/3 if to be transformed to a 5-poit scale
                                      fltlnl %in% c(7,8,9) ~ NA),
         id = idno,
         age = ifelse(agea == 999, NA, agea),
         edu_years = eduyrs,
         education = factor(case_when(eisced %in% c(55,77,88, 99) ~ NA, # ISCED-based
                                      eisced == 1 ~ "primary",
                                      eisced %in% c(2,3,4) ~ "secondary",
                                      eisced %in% c(5,6) ~ "post-secondary vocational or lower tertiary",
                                      eisced == 7 ~ "upper tertiary"), # master's degree or higher
                            levels = c("primary", "secondary", "post-secondary vocational or lower tertiary", "upper tertiary"), 
                            ordered = T), 
         urbanicity = case_when(domicil %in% c(1,2,3) ~ "urban", # a town or a city
                                domicil %in% c(4,5) ~ "rural"),
         household_composition = case_when(hhmmb == 1 ~ "living alone",
                                           hhmmb %in% c(77, 88, 99) ~ NA,
                                           hhmmb > 1 & hhmmb < 50 ~ "other"),
         marital_status = case_when(maritalb %in% c(1,2) ~ "married/cohabitating/registered union",
                                    maritalb %in% c(3,4,5) ~ "separated/divorced/widowed",
                                    maritalb == 6 ~ "never married",
                                    maritalb %in% c(77,88,99) ~ NA),
         work_status = factor(case_when(mnactic %in% c(1,7) ~ "employed",
                                        mnactic %in% c(3,4) ~ "unemployed",
                                        mnactic %in% c(2,5,6,8) ~ "economically inactive", 
                                        mnactic %in% c(9,77,88, 99) ~ NA), 
                              levels = c("employed", "unemployed", "economically inactive"), ordered = T),
         household_income = as.numeric(case_when(hinctnta %in% c(77,88,99) ~ NA, # in country-specific deciles, from 1 to 10
                                                 TRUE ~ hinctnta)),
         gender = case_when(gndr == 1 ~ "male", 
                            gndr == 2 ~ "female"),
         country = case_when(cntry == "AT" ~ "Austria",
                             cntry == "BE" ~ "Belgium",
                             cntry == "BG" ~ "Bulgaria",
                             cntry == "CH" ~ "Switzerland",
                             cntry == "CY" ~ "Cyprus",
                             cntry == "DE" ~ "Germany",
                             cntry == "EE" ~ "Estonia",
                             cntry == "ES" ~ "Spain",
                             cntry == "FI" ~ "Finland",
                             cntry == "FR" ~ "France",
                             cntry == "GB" ~ "United Kingdom",
                             cntry == "GR" ~ "Greece",
                             cntry == "HR" ~ "Croatia",
                             cntry == "HU" ~ "Hungary",
                             cntry == "IE" ~ "Ireland",
                             cntry == "IS" ~ "Iceland",
                             cntry == "IL" ~ "Israel",
                             cntry == "IT" ~ "Italy",
                             cntry == "LT" ~ "Lithuania",
                             cntry == "LV" ~ "Latvia",
                             cntry == "ME" ~ "Montenegro",
                             cntry == "NL" ~ "Netherlands",
                             cntry == "NO" ~ "Norway",
                             cntry == "PL" ~ "Poland",
                             cntry == "PT" ~ "Portugal",
                             cntry == "RS" ~ "Serbia",
                             cntry == "SE" ~ "Sweden",
                             cntry == "SI" ~ "Slovenia",
                             cntry == "SK" ~ "Slovakia"),
         smoking = case_when(cgtsmok %in% c(1,2,3) ~ "daily/occasional",
                             cgtsmok %in% c(4,5,6) ~ "none",
                             cgtsmok %in% c(7,8,9) ~ NA),
         body_weight = ifelse(weighta %in% c(777,888,999), NA, weighta),
         height = ifelse(height %in% c(777,888,999), NA, height)/100,
         BMI = body_weight/height^2,
         obesity = ifelse(BMI >= 30, 1, 0),
         phys_activity = case_when(dosprt %in% c(77,88,99) ~ NA,
                                   TRUE ~ dosprt),
         weight = pspwght, # a post-stratification weight
         psu, # primary sampling unit - clustering
         stratum, # stratification
         dataset = "ESS11") %>%  
  select(dataset, country, id, weight, loneliness, conti_dir_loneli,
         age, gender, education, household_composition, marital_status, work_status, household_income, urbanicity,
         smoking, BMI, obesity, phys_activity,
         psu, stratum)


# European Social Survey 7 (ESS7, 2014, probability-based household sample, in-person interviewed)
ESS7_data <- ESS7 %>% 
  mutate(loneliness = case_when(fltlnl %in% c(3,4) ~ 1, # all or most of time
                                fltlnl %in% c(1,2) ~ 0, # none or some of time
                                fltlnl %in% c(7,8,9) ~ NA), # refusal, don't know or no answer
         conti_dir_loneli = case_when(fltlnl %in% c(1,2,3,4) ~ fltlnl, # ~ 1 + (fltlnl - 1) * 4/3
                                      fltlnl %in% c(7,8,9) ~ NA),
         id = idno,
         country = case_when(cntry == "AT" ~ "Austria",
                             cntry == "BE" ~ "Belgium",
                             cntry == "CH" ~ "Switzerland",
                             cntry == "CZ" ~ "Czechia",
                             cntry == "DE" ~ "Germany",
                             cntry == "DK" ~ "Denmark",
                             cntry == "EE" ~ "Estonia",
                             cntry == "ES" ~ "Spain",
                             cntry == "FI" ~ "Finland",
                             cntry == "FR" ~ "France",
                             cntry == "GB" ~ "United Kingdom",
                             cntry == "HR" ~ "Croatia",
                             cntry == "HU" ~ "Hungary",
                             cntry == "IE" ~ "Ireland",
                             cntry == "IS" ~ "Iceland",
                             cntry == "IL" ~ "Israel",
                             cntry == "IT" ~ "Italy",
                             cntry == "LT" ~ "Lithuania",
                             cntry == "LV" ~ "Latvia",
                             cntry == "ME" ~ "Montenegro",
                             cntry == "NL" ~ "Netherlands",
                             cntry == "NO" ~ "Norway",
                             cntry == "PL" ~ "Poland",
                             cntry == "PT" ~ "Portugal",
                             cntry == "RS" ~ "Serbia",
                             cntry == "SE" ~ "Sweden",
                             cntry == "SI" ~ "Slovenia",
                             cntry == "SK" ~ "Slovakia"),
         age = ifelse(agea == 999, NA, agea),
         edu_years = eduyrs,
         education = factor(case_when(eisced %in% c(55,77,88, 99) ~ NA, # ISCED-based
                                      eisced == 1 ~ "primary",
                                      eisced %in% c(2,3,4) ~ "secondary",
                                      eisced %in% c(5,6) ~ "post-secondary vocational or lower tertiary",
                                      eisced == 7 ~ "upper tertiary"), # master's degree or higher
                            levels = c("primary", "secondary", "post-secondary vocational or lower tertiary", "upper tertiary"), 
                            ordered = T), 
         household_composition = case_when(hhmmb == 1 ~ "living alone",
                                           hhmmb %in% c(77, 88, 99) ~ NA,
                                           hhmmb > 1 & hhmmb < 50 ~ "other"),
         marital_status = case_when(maritalb %in% c(1,2) ~ "married/cohabitating/registered union",
                                    maritalb %in% c(3,4,5) ~ "separated/divorced/widowed",
                                    maritalb == 6 ~ "never married",
                                    maritalb %in% c(77,88,99) ~ NA),
         urbanicity = case_when(domicil %in% c(1,2,3) ~ "urban", # a town or a city
                                domicil %in% c(4,5) ~ "rural"),
         work_status = factor(case_when(mnactic %in% c(1,7) ~ "employed",
                                        mnactic %in% c(3,4) ~ "unemployed",
                                        mnactic %in% c(2,5,6,8) ~ "economically inactive", 
                                        mnactic %in% c(9,77,88, 99) ~ NA),
                              levels = c("employed", "unemployed", "economically inactive"), ordered = T),
         household_income = as.numeric(case_when(hinctnta %in% c(77,88,99) ~ NA, # in country-specific deciles, from 1 to 10
                                                 TRUE ~ hinctnta)),
         gender = case_when(gndr == 1 ~ "male", 
                            gndr == 2 ~ "female"),
         smoking = case_when(cgtsmke %in% c(1,2) ~ "daily/occasional",
                             cgtsmke %in% c(3,4,5) ~ "none",
                             cgtsmke %in% c(7,8,9) ~ NA),
         body_weight = ifelse(weight %in% c(777,888,999), NA, weight),
         height = ifelse(height %in% c(777,888,999), NA, height)/100,
         BMI = body_weight/height^2,
         obesity = ifelse(BMI >= 30, 1, 0),
         phys_activity = case_when(dosprt %in% c(77,88,99) ~ NA,
                                   TRUE ~ dosprt),
         weight = pspwght, # a post-stratification weight
         psu, # primary sampling unit - clustering
         stratum, # stratification
         dataset = "ESS7") %>%  
  select(dataset, country, id, weight, loneliness, conti_dir_loneli,
         age, gender, education, household_composition, marital_status, work_status, household_income, urbanicity,
         smoking, BMI, obesity, phys_activity,
         psu, stratum)


# EU Loneliness Survey: EU27 sample (2022, agency online panel)
EU27_data <- EU27 %>%                                        # In general, how lonely do you feel?
  mutate(loneliness = case_when(loneliness_direct %in% c(1,2) ~ 1, # all or most of time
                                loneliness_direct %in% c(3,4,5) ~ 0, # some, a little, or none of time
                                loneliness_direct %in% c(998,999) ~ NA), # don't know or prefer not to say
         #conti_dir_loneli = case_when(loneliness_direct %in% c(1,2,3,4,5) ~ 6 - loneliness_direct,
         #                             loneliness_direct %in% c(998,999) ~ NA), # 5-point scale
         conti_dir_loneli = case_when(loneliness_direct == 1 ~ 4,
                                      loneliness_direct == 2 ~ 3,
                                      loneliness_direct == 3 ~ 2,
                                      loneliness_direct %in% c(4,5) ~ 1,
                                      loneliness_direct %in% c(998,999) ~ NA),
         across(starts_with("loneliness_ucla"), ~ ifelse(.x == 999, NA, .x)), # prefer not to say
         across(starts_with("loneliness_djg"), ~ ifelse(.x == 999, NA, .x)),
         across(c(loneliness_djg_a:loneliness_djg_c), ~ 4 - .x)) %>%  
  mutate(ucla_loneli = rowMeans(select(., loneliness_ucla_a:loneliness_ucla_c), na.rm = T),
         djg_loneli = rowMeans(select(., loneliness_djg_a:loneliness_djg_f), na.rm = T),
         across(c(ucla_loneli, djg_loneli), ~ (.x - 1)/(3 - 1) * 100),
         id = id_pers, 
         country = as.character(as_factor(country)),
         dataset = "EU27",
         gender = case_when(gender == 1 ~ "male", 
                            gender == 2 ~ "female"),
         age = ifelse(age == 997, NA, age),
         education = factor(case_when(education == 999 ~ NA,
                                      education %in% c(1,2) ~ "primary",
                                      education == 3 ~ "secondary",
                                      education == 4 ~ "post-secondary vocational or lower tertiary",
                                      education == 5 ~ "upper tertiary"),
                            levels = c("primary", "secondary", "post-secondary vocational or lower tertiary", "upper tertiary"), 
                            ordered = T), 
         household_composition = case_when(adults_n == 1 & children_hh == 0 & others_n == 0 ~ "living alone",
                                           !(adults_n == 1 & children_hh == 0 & others_n == 0) ~ "other"),
         marital_status = case_when(relationship %in% c(1,2) ~ "never married",
                                    relationship == 3 ~ "married/cohabitating/registered union",
                                    relationship %in% c(4,5) ~ "separated/divorced/widowed",
                                    relationship == 999 ~ NA),
         work_status = factor(case_when(work_status == 1 ~ "employed",
                                        work_status %in% c(3,4) ~ "unemployed",
                                        work_status %in% c(2,5,6,7) ~ "economically inactive", 
                                        work_status == 8 ~ NA),
                              levels = c("employed", "unemployed", "economically inactive"), ordered = T),
         household_income = as.numeric(case_when(income_decile == 999 ~ NA, # in country-specific deciles, from 1 to 10
                                                 TRUE ~ income_decile)),
         urbanicity = case_when(municipality_type == 1 ~ "rural",
                                municipality_type %in% c(2,3) ~ "urban"),
         smoking = case_when(smoking %in%  c(1,2) ~ "daily/occasional",
                             smoking == 3 ~ 'none',
                             smoking == 999 ~ NA),
         body_weight = ifelse(country == "Ireland" & weight__open_eire > 110, weight__open_eire*0.45, weight__open), # pounds to kilos
         height = ifelse(country == "Ireland" & height__open_eire < 120, height__open_eire*2.54/100, height__open/100), # inches to cm
         BMI = body_weight/height^2,
         obesity = ifelse(BMI >= 30, 1, 0),
         phys_activity = case_when(exercise %in% c(998, 999) ~ NA,
                                   TRUE ~ exercise - 1),
         weight = w_country) %>% 
  select(dataset, country, id, weight, loneliness, conti_dir_loneli, ucla_loneli, djg_loneli,
         age, gender, education, household_composition, marital_status, work_status, household_income, urbanicity,
         smoking, BMI, obesity, phys_activity)


# EU Loneliness Survey: EU4 sample (probability-based, online interviewing)
EU4_data <- EU4 %>%                                          # How much of the time, during the past 4 weeks, have you been feeling lonely?
  mutate(loneliness = case_when(loneliness_direct %in% c(1,2) ~ 1, # all or most of time
                                loneliness_direct %in% c(3,4,5) ~ 0, # some, a little, or none of time
                                loneliness_direct %in% c(998,999) ~ NA), # don't know or prefer not to say
         conti_dir_loneli = case_when(loneliness_direct == 1 ~ 4,
                                      loneliness_direct == 2 ~ 3,
                                      loneliness_direct == 3 ~ 2,
                                      loneliness_direct %in% c(4,5) ~ 1,
                                      loneliness_direct %in% c(998,999) ~ NA),
         across(starts_with("loneliness_ucla"), ~ ifelse(.x == 999, NA, .x)), # prefer not to say
         across(starts_with("loneliness_djg"), ~ ifelse(.x == 999, NA, .x)),
         across(c(loneliness_djg_a:loneliness_djg_c), ~ 4 - .x)) %>%  
  mutate(ucla_loneli = rowMeans(select(., loneliness_ucla_a:loneliness_ucla_c), na.rm = T),
         djg_loneli = rowMeans(select(., loneliness_djg_a:loneliness_djg_f), na.rm = T),
         across(c(ucla_loneli, djg_loneli), ~ (.x - 1)/(3 - 1) * 100),
         id = id, 
         gender = case_when(gender == 1 ~ "male", 
                            gender == 2 ~ "female"),
         age = age,
         education = factor(case_when(education == 999 ~ NA,
                                      education %in% c(1,2) ~ "primary",
                                      education %in% c(3,4) ~ "secondary",
                                      education %in% c(5,6,7) ~ "post-secondary vocational or lower tertiary",
                                      education %in% c(8,9) ~ "upper tertiary"),
                            levels = c("primary", "secondary", "post-secondary vocational or lower tertiary", "upper tertiary"), 
                            ordered = T), 
         household_composition = case_when(adults_n == 1 & children_hh == 0 & others_n == 0 ~ "living alone",
                                           !(adults_n == 1 & children_hh == 0 & others_n == 0) ~ "other"),
         marital_status = case_when(relationship %in% c(1,2) ~ "never married",
                                    relationship == 3 ~ "married/cohabitating/registered union",
                                    relationship %in% c(4,5) ~ "separated/divorced/widowed",
                                    relationship == 999 ~ NA),
         work_status = factor(case_when(work_status == 1 ~ "employed",
                                        work_status %in% c(3,4) ~ "unemployed",
                                        work_status %in% c(2,5,6,7) ~ "economically inactive", 
                                        work_status == 8 ~ NA),
                              levels = c("employed", "unemployed", "economically inactive"), ordered = T),
         household_income = as.numeric(case_when(income_month %in% c(998,999) ~ NA, # in country-specific deciles, from 1 to 10
                                                 TRUE ~ income_month)),
         urbanicity = case_when(municipality_type == 1 ~ "rural",
                                municipality_type %in% c(2,3) ~ "urban"),
         smoking = case_when(smoking %in%  c(1,2) ~ "daily/occasional",
                             smoking == 3 ~ "none",
                             smoking == 999 ~ NA),
         phys_activity = case_when(exercise %in% c(998, 999) ~ NA,
                                   TRUE ~ exercise - 1),
         weight = w_country,
         country = as.character(as_factor(country)),
         dataset = "EU4") %>% # country-specific weights
  select(dataset, country, id, weight, loneliness, conti_dir_loneli, ucla_loneli, djg_loneli,
         age, gender, education, household_composition, marital_status, work_status, household_income, urbanicity,
         smoking, phys_activity)




##############################
## DESCRIPTIVE STATISTICS ###
#############################

ESS_EULS_dataset <- bind_rows(ESS7_data, ESS11_data, EU27_data, EU4_data) %>% 
  mutate(across(c(gender, urbanicity, household_composition, marital_status, smoking, obesity, loneliness),
                ~ as.factor(.x)),
         across(c(education, work_status),
                ~ factor(.x, ordered = F)),
         across(c(age, household_income, BMI, phys_activity, ucla_loneli, djg_loneli, conti_dir_loneli),
                ~ as.numeric(.x)))
write.csv(ESS_EULS_dataset, "C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Datasets/Data/ESS_EULS_dataset.csv")


# SAMPLE SIZES
ESS_EULS_sample_sizes <- ESS_EULS_dataset  %>% 
  select(dataset, country, id) %>% 
  group_by(dataset, country) %>% 
  count(dataset, country) %>% 
  pivot_wider(names_from = dataset, values_from = n) %>% 
  select(country, ESS7, ESS11, EU27, EU4) %>% 
  arrange(country)
write.csv(ESS_EULS_sample_sizes, "C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Prevalence/Results/Tables/ESS_EULS_sample_sizes.csv")


# SAMPLE CHARACTERISTICS
ESS_EULS_props <- ESS_EULS_dataset  %>% 
  dplyr::select(dataset, country, id, weight, gender, education, household_composition, marital_status, work_status, urbanicity, smoking) %>% 
  pivot_longer(-c(dataset, country, id, weight)) %>% 
  group_by(dataset, country, name) %>%
  summarise(freq = Hmisc::wtd.table(x = value, weights = weight, na.rm = T, type = "table"),
            prop = weights::wpct(x = value, weight = weight, na.rm = T),
            stat = paste(formatC(round(freq, 1), format = "f", digits = 1), 
                         paste0("(", formatC(round(prop*100, 1), format = "f", digits = 1),  ")"))) %>% 
  ungroup() %>% 
  mutate(value = names(freq)) %>% 
  dplyr::select(dataset, country, name, value, stat) %>% 
  pivot_wider(names_from = dataset, values_from = stat) %>% 
  select(country, name, value, ESS7, ESS11, EU27, EU4) 

ESS_EULS_means <- ESS_EULS_dataset  %>% 
  dplyr::select(dataset, country, id, weight, age, household_income, phys_activity, BMI) %>% 
  pivot_longer(-c(dataset, country, id, weight)) %>% 
  group_by(dataset, country, name) %>%
  summarise(mean = formatC(round(weighted.mean(value, w = weight, na.rm = T), 1), format = "f", digits = 1),
            sd = formatC(round((Hmisc::wtd.var(value, weights = weight, na.rm = T)/sqrt(Hmisc::wtd.var(value, weights = weight, na.rm = T))), 1), format = "f", digits = 1)) %>% 
  mutate(stat = paste0(mean, " (", sd, ")")) %>% 
  select(dataset, country, name, stat) %>% 
  pivot_wider(names_from = dataset, values_from = stat) %>% 
  mutate(value = NA) %>% 
  select(country, name, value, ESS7, ESS11, EU27, EU4) 

ESS_EULS_sample_characteristics <- rbind(ESS_EULS_props, ESS_EULS_means) %>%
  ungroup() %>% 
  mutate(name = factor(name, ordered = T,
                       levels = c("age", "gender", "education", "household_composition", "marital_status", "work_status", "household_income", "urbanicity",
                                  "phys_activity", "smoking", "BMI"))) %>% 
  filter(!value %in% c("rural", "none", "other")) %>% 
  arrange(name, value, country)
write.csv(ESS_EULS_sample_characteristics, "C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Prevalence/Results/Tables/ESS_EULS_sample_characteristics.csv")


overall_props <- ESS_EULS_dataset  %>% 
  dplyr::select(id, weight, gender, household_composition, marital_status, education, work_status, urbanicity, smoking) %>% 
  pivot_longer(-c(id, weight)) %>% 
  group_by(name) %>%
  summarise(freq = Hmisc::wtd.table(x = value, weights = weight, na.rm = T, type = "table"),
            prop = weights::wpct(x = value, weight = weight, na.rm = T),
            stat = paste(formatC(round(freq, 1), format = "f", digits = 1), 
                         paste0("(", formatC(round(prop*100, 1), format = "f", digits = 1),  ")"))) %>% 
  ungroup() %>% 
  mutate(value = names(freq)) %>% 
  dplyr::select(name, value, stat) 

overall_means <- ESS_EULS_dataset  %>% 
  dplyr::select(id, weight, age, household_income, phys_activity, BMI) %>% 
  pivot_longer(-c(id, weight)) %>% 
  group_by(name) %>%
  summarise(mean = formatC(round(weighted.mean(value, w = weight, na.rm = T), 1), format = "f", digits = 1),
            sd = formatC(round((Hmisc::wtd.var(value, weights = weight, na.rm = T)/sqrt(Hmisc::wtd.var(value, weights = weight, na.rm = T))), 1), format = "f", digits = 1)) %>% 
  mutate(stat = paste0(mean, " (", sd, ")")) %>% 
  select(name, stat) %>% 
  mutate(value = NA) %>% 
  select(name, value, stat) 

overall_sample_characteristics <- rbind(overall_props, overall_means) %>%
  ungroup() %>% 
  mutate(name = factor(name, ordered = T,
                       levels = c("age", "gender", "education", "household_composition", "marital_status", "work_status", "household_income", "urbanicity",
                                  "phys_activity", "smoking", "BMI"))) %>% 
  filter(!value %in% c("rural", "none", "other")) %>% 
  arrange(name, value)
write.csv(ESS_EULS_sample_characteristics, "C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Prevalence/Results/Tables/ESS_EULS_sample_characteristics.csv")

# age range of the study samples: from 16 to 80 years
ESS_EULS_age_ranges <- ESS_EULS_dataset  %>% 
  dplyr::select(dataset, country, id, age) %>% 
  group_by(dataset, country) %>%
  summarise(min = min(age, na.rm = T),
            max = max(age, na.rm = T)) %>% 
  mutate(stat = paste0(min, " to ", max)) %>% 
  select(dataset, country, stat) %>% 
  pivot_wider(names_from = dataset, values_from = stat) %>% 
  select(country, ESS7, ESS11, EU27, EU4) %>% 
  arrange(country)
write.csv(ESS_EULS_age_ranges, "C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Prevalence/Results/Tables/ESS_EULS_age_ranges.csv")

# MISSING DATA
ESS_EULS_missing_prop_lonely <- ESS_EULS_dataset %>% 
  dplyr::select(dataset, country, id, weight, loneliness) %>% 
  mutate(loneliness = ifelse(is.na(loneliness), 1, 0)) %>%   # missing (NA) recoded into 1
  group_by(dataset, country) %>%
  summarise(freq = Hmisc::wtd.table(x = loneliness, weights = weight, na.rm = F, type = "table"),
            prop = weights::wpct(x = loneliness, weight = weight, na.rm = F)*100,
            stat = formatC(round(prop, 1), format = "f", digits = 1)) %>% 
  ungroup() %>% 
  mutate(value = names(freq)) %>% 
  filter(value == "1") %>% 
  dplyr::select(dataset, country, stat) %>% 
  pivot_wider(names_from = dataset, values_from = stat) %>% 
  select(country, ESS7, ESS11, EU27, EU4) %>% 
  arrange(country)
write.csv(ESS_EULS_missing_prop_lonely, "C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Prevalence/Results/Tables/ESS_EULS_missing_prop_lonely.csv")




############################
### MULTIPLE IMPUTATION ###
###########################

ESS_EULS_dataset <- ESS_EULS_dataset %>% 
  mutate(education = factor(education, ordered = T),
         -c(ucla_loneli, djg_loneli, conti_dir_loneli))

init_imp <- mice(ESS_EULS_dataset, maxit = 0) # dry run for specifying predictors and imputation model

pred_matrix <- init_imp$predictorMatrix
pred_matrix[c("id","weight"),] <- 0 # not imputing IDs and post-stratification weights
pred_matrix[, c("id","weight")] <- 0

imp_method <- init_imp$method
# predictive mean matching (pmm) for continuous,
# ordinal/proportional odds logistic (polr) for ordered categorical,
# logistic regression (logreg) for binary,
# polytomous regression (polyreg) for categorical

set.seed(123321)
imp_data <- mice(ESS_EULS_dataset,
                 m = 5, # number of imputed datasets
                 maxit = 10, # number of iterations
                 method = imp_method,
                 predictorMatrix = pred_matrix)

plot(imp_data)
imp_data_long <- complete(imp_data, "long", include = T)
write.csv(imp_data_long, "C:/Users/hp/Desktop/PhD_Oslo/Research/Loneliness_Norms_Contexts/Datasets/Data/imp_data_long.csv")
