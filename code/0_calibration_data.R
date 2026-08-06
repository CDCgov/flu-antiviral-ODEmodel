## Calibration data -----------

# Disease burden estimates, from https://www.cdc.gov/flu-burden/php/data-vis/past-seasons.html
# Accessed 5 Feb 2026
db_tmp <- read.csv(here(filepath, 'data/flu-disease-burden-past-season-estimates.csv')) %>%
  clean_names() 

# All ages
db <- db_tmp %>% filter(age_group == "All Ages") %>%
  select(season = flu_season, contains("symptomatic"), contains("hospitalizations")) 


# Add season severity classifications, from https://www.cdc.gov/flu/php/surveillance/past-seasons.html
# Accessed 28 Dec 2025
db$severity <- c("moderate", "low", "moderate", "moderate", 
                 "high", "moderate", "moderate", "high", "moderate", 
                 "moderate", "low", "moderate", "moderate", "high")


# By age group 
db_age <- db_tmp %>% filter(age_group != "All Ages") %>%
  select(season = flu_season, age_txt = age_group, contains("symptomatic"), contains("hospitalizations")) %>%
  left_join(db %>% distinct(season, severity)) %>%
  mutate(age_txt = factor(age_txt, levels = c("0-4 yrs", "5-17 yrs", "18-49 yrs", "50-64 yrs", "65+ yrs")))


