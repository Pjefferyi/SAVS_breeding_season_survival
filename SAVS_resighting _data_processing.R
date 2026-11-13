# This script is to process the resighting data collected during the breeding season 

# Load libraries 
library(tidyverse)
library(sf)
library(ggplot2)
library(anytime)
library(patchwork)

# Path to Jelany's local data ---- 
local.data <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/PhD_Thesis/Data"

# Path to folder for figures ----
fig.dir <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/PhD_Thesis/Chapter_3_breeding_survival/breeding_survival_figures"

# Load all required data ----

## Load the resighting data from breeding surveys ---- 
br.survey <- read_csv(paste0(local.data,"/breeding_survival/SAVS_resightings.csv")) %>%
  mutate(date = anytime(date)) %>% filter(likely_band != "unknown") %>%
  group_by(likely_band) %>%
  mutate(first_sighting = last(date)) %>%
  ungroup() %>%
  mutate(likely_band = fct_reorder(factor(likely_band), as.numeric(rev(first_sighting))))

## Load morphological data ----
morph <- read_csv(paste0(local.data,"/Kent_morpho_nest_ind/SAVS_dataset_backup_20260929/JD_morphology_20260929.csv"),
                  col_types = cols(
                    band_id = col_character(),
                    location = col_character(),
                    weight = col_double(),
                    tarsus = col_double(),
                    wing_l = col_double(),
                    tail_l = col_double(),
                    bill_l = col_double(),
                    bill_d = col_double(),
                    date_measured = col_date(),
                    time_measured = col_double(),
                    sex = col_character(),
                    age_in_yrs = col_double(),
                    age_in_d = col_double(),
                    color = col_character(),
                    year = col_double(),
                    fat_score = col_double(),
                    bander_id = col_character(),
                    Notes = col_character()
                  )) %>% mutate(date_measured = anytime(ymd(date_measured)))

## Load Individual data ---- 
individs <- read_csv(paste0(local.data,"/Kent_morpho_nest_ind/SAVS_dataset_backup_20260929/JD_individuals_20260929.csv"),
                     col_types = cols(
                       band_id = col_character(),
                       adult_color = col_character(),
                       hy_color = col_character(),
                       sex = col_character(),
                       birth_year = col_double(),
                       known_age = col_logical(),
                       banding_age = col_character(),
                       mother_id = col_character(),
                       father_id = col_character(),
                       banding_date = col_date(),
                       natal_nestkey = col_character(),
                       mother_color = col_character(),
                       father_color = col_character(),
                       Notes = col_character()
                     ))%>% mutate(banding_date = anytime(ymd(banding_date)))

## Load Nest data ----
nests <- read_csv(paste0(local.data,"/Kent_morpho_nest_ind/SAVS_dataset_backup_20260929/JD_nests_20260929.csv"),
                  col_types = cols(
                    nest_key = col_character(),
                    year = col_double(),
                    nest_id = col_character(),
                    female_id = col_character(),
                    male_id = col_character(),
                    eggs_laid = col_double(),
                    found_date = col_date(),
                    hatch_date = col_date(),
                    eggs_hatched = col_double(),
                    fledglings = col_double(),
                    fledge_date = col_date(),
                    study_site = col_character(),
                    fate = col_character(),
                    mating_status = col_character(),
                    female_number = col_character(),
                    female_color = col_character(),
                    male_color = col_character(),
                    clutch_number = col_character(),
                    x_coord = col_double(),
                    y_coord = col_double(),
                    fail_date = col_character(),
                    experiment = col_character(),
                    Notes = col_character()
                  ))


# Add fall banding recaptures to the resightings ----

# filter the morphology dataaset to get capture data
fall.captures <- morph %>% 
  filter(year > 2024) %>%
  mutate(j_date_measured = yday(date_measured)) %>%
  filter(j_date_measured >= 227 & j_date_measured < 258) %>%
  group_by(band_id) %>%
  filter(date_measured == first(date_measured)) %>%
  mutate(coll_date = update(date_measured, month = 9, day = 1)) 

# We will format these capture to match the breeding survey
fall.survey <-fall.captures %>% 
  filter(color %in%  br.survey$likely_band) %>%
  select(band_id, color, coll_date) %>%
  rename(likely_band = color, date = coll_date) %>%
  mutate(period = "postbreeding", site = "S")

# merge fall captures with breeding dataset
br.survey <- bind_rows(br.survey, fall.survey)

# Initial plot of the re sightings ----
ggplot() + 
  geom_tile(data = br.survey, aes(x= as.factor(date), y = likely_band, fill = period), colour = "white") + 
  labs(fill = "Period") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 15))

# We add individual data to the resightings dataset ----
individuals <- read_csv(paste0(local.data,"/Kent_morpho_nest_ind/SAVS_dataset_backup_20260929/JD_individuals_20260929.csv"),
                  col_types = cols(
                    band_id = col_character(),
                    adult_color = col_character(),
                    hy_color = col_character(),
                    sex = col_character(),
                    birth_year = col_double(),
                    known_age = col_logical(),
                    banding_age = col_character(),
                    mother_id = col_character(),
                    father_id = col_character(),
                    banding_date = col_date(),
                    natal_nestkey = col_character(),
                    mother_color = col_character(),
                    father_color = col_character(),
                    Notes = col_character()
                  ))%>% mutate(banding_date = anytime(banding_date)) %>%
  group_by(adult_color) %>%
  # Here we assume that the birds that have been are the most recent birds to be using their specific combos 
  # This would break down if any combos have been re-used while a bird is still alive, or if this analysis is performed again in the future
  # So check the data 
  filter(banding_date == max(banding_date))

br.survey.ex <- br.survey %>% select(-band_id) %>% merge(individuals, by.x = "likely_band", by.y = "adult_color") %>%
  arrange(desc(banding_date)) %>% mutate(likely_band = factor(likely_band, levels = unique(likely_band), ordered = T)) %>% filter(sex %in% c("F", "M"))

# Re-sightings by period ----
ggplot() + 
  geom_tile(data = br.survey.ex, aes(x= as.factor(date), y = likely_band, fill = period), colour = "white")+ 
  labs(fill = "Period") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 15)) +
  facet_wrap(~sex)

View(br.survey.ex)

# Define periods by bird ---

# Parameters: 
# prep time: estimated at 17 days - amount of days necessary to build a nest, lay eggs, and incubate
# If a  bird's nest failed without hatching, we estimate that the females was a the halfway point of the process at the tiem of the failure, i.e., 7 days
# nestling care: fledge date - hatch date (generally about 11 days) 
# fledgling care: estimated at 13 days 

# We create a dataframe where each row is a day in the breeding season of a bird
# A bird is considered to be in the pre-breeding period until the 17 days before it's first nest
# A bird enters the preparation phase 17 days before the eggs hatch
# A birds enters the nestling care period until the young fledge
# A bird then enters the fledgling care period until: a) 13 days have passed or b) it enters the preparation phase for a new nest
# Once the fledgling care period ends for the last nest of a bird, it enters the postbreedign period. 

# We will create a dataframe with bird periods

# We begin by adding the nestling and fledgling care periods 
# We then add the preparation period
# Lastly, we define the pre-breeding and postbreeding period before and after the first and last nest
# any gaps that remain will be preparation

# first we will modify the nest database to calculate some breeding periods ----
nests.br.periods <- nests %>% filter(year >= 2025) %>% group_by(nest_key) %>%
  mutate(nestling_start = case_when(fate == "FL" ~ hatch_date,
                                   fate == "AN" ~ hatch_date,
                                   fate == "AE" ~ NA,
                                   fate == "U" ~ NA),
         nestling_end = case_when(fate == "FL" ~ fledge_date,
                                         fate == "AN" ~ anytime(fail_date),
                                         fate == "AE" ~ NA,
                                         fate == "U" ~ NA),
         prep_start = case_when(fate == "FL" ~ anytime(hatch_date) - days(17),
                                     fate == "AN" ~ anytime(hatch_date) - days(17),
                                     fate == "AE" ~ anytime(fail_date) - days(8),
                                     fate == "U" ~ NA),
         prep_end = case_when(fate == "FL" ~ anytime(hatch_date),
                                     fate == "AN" ~ anytime(hatch_date),
                                     fate == "AE" ~ anytime(fail_date),
                                     fate == "U" ~ NA),
         fledgling_start =  case_when(fate == "FL" ~ anytime(fledge_date)),
         fledgling_end =  case_when(fate == "FL" ~ anytime(fledge_date) + days(13))) %>%
  # next we pivot this 
  pivot_longer(
    cols = c(
      nestling_start, nestling_end, prep_start,
      fledgling_start, fledgling_end, prep_end,
    ),
    names_to = c("care_type", ".value"),
    names_pattern = "(nestling|fledgling|prep)_(start|end)"
  ) %>%
  rename(
    start_time = start,
    end_time = end
  ) %>%
  mutate(
    state = recode(
      care_type,
      nestling = "nestling care",
      fledgling = "fledgling care",
      prep = "preparation",
    )
  ) %>%
  filter(!is.na(start_time) & !is.na(end_time) & ! is.na(female_id)) %>%
  select(nest_key, male_id, male_color,female_id, female_color, state, start_time, end_time, year)

# Save breeding state data ----
write_csv(nests.br.periods, paste0(local.data,"/breeding_survival/nesting_states.csv"))

# We can apply this state information to the resighting survey data ----

# We estimate the timing of the pre-breeding state
# The length of the post0breeding state is not recoverable because we do not know when and if individuals departed Kent Island for fall migration 

# We will expand the dataset (with complete())and assign states even to dates when no survey attempts where made 
# Then we can calculate the amount of time spent in each stage. 

# Here we start with females 
br.survey.states.f <- br.survey.ex %>%
  rename("survey_start_time" = "start_time", 
         "survey_end_time" = "end_time") %>%
  filter(sex == "F") %>%
  mutate(year = year(date)) %>%
  # we expand the resighting dataset to include even days when surveys were not conducted. 
  # This is to ensure that we can make daily state assignments
  group_by(band_id, year) %>%
  mutate(date = as.Date(date)) %>%
  complete(date = seq.Date(min(date), max(date), by = "day"), nesting(likely_band), fill = list(period = "Not seen or surveyed")) %>%
  group_by(date)%>%
  fill(survey_start_time, survey_end_time, surveyor, precip, cloud_cover, wind, fog, Site, .direction = "down") %>%
  ungroup() %>%
  # add 1s and 2s to indicate whether an individual was resighted 
  mutate(sighted = ifelse(is.na(band_as_seen), 1, 2)) %>%
  # merge with the nesting peruiod data
  left_join(
    nests.br.periods,
    by = join_by("band_id" == "female_id", "year" == "year"),
                 relationship = "many-to-many") %>%
  # Check if each female is in any of the breeding states 
  mutate(in_state = !is.na(start_time) &
      !is.na(end_time) &
      date >= start_time &
      date < end_time) %>%
  # If a bird is in none of the states, we will give it the "preparation" state
  group_by(band_id, date) %>%
  filter(in_state == T | all(in_state == F)) %>%
  mutate(state = ifelse(in_state == F, "preparation", state)) %>%
  distinct(band_id, date, .keep_all = T) %>%
 # then we can define pre-breeding and post-breeding states 
  group_by(band_id, year) %>%
  mutate(breeding.start = min(start_time),
         breeding.end = max(end_time),
         state = case_when(!is.na(breeding.start) & date <= breeding.start ~ "prebreeding",
                           !is.na(breeding.end) & date > breeding.end ~ "postbreeding",
                           is.na(breeding.start) & yday(date) <= 193 ~ "prebreeding",
                           is.na(breeding.end) & yday(date) > 227 ~ "postbreeding",
                           is.na(breeding.end) & yday(date) > 227 & yday(date) < 227 ~ "preparation",
                           .default = state)) %>%
  arrange(band_id, date) %>%
  mutate(period = factor(period, levels = c("prebreeding", "breeding", "postbreeding")),
         state = factor(state, levels = c("prebreeding", "preparation", "nestling care", "fledgling care", "postbreeding"))) %>%
  # lastly, we will modify the "start-time" and "end_time" column so they include updated information 
  group_by(band_id) %>%
  mutate(state_id = cumsum(state != lag(state, default = first(state)))) %>%
  group_by(band_id, state_id) %>%
  mutate(start_time = min(date),
         end_time = max(date)) %>%
  ungroup()

# Now we can assign states to males
br.survey.states.m <- br.survey.ex %>%
  rename("survey_start_time" = "start_time", 
         "survey_end_time" = "end_time") %>%
  filter(sex == "M") %>%
  mutate(year = year(date)) %>%
  # we expand the resighting dataset to include even days when surveys were not conducted. 
  # This is to ensure that we can make daily state assignments
  group_by(band_id, year) %>%
  mutate(date = as.Date(date)) %>%
  complete(date = seq.Date(min(date), max(date), by = "day"), nesting(likely_band), fill = list(period = "Not seen or surveyed")) %>%
  group_by(date)%>%
  fill(survey_start_time, survey_end_time, surveyor, precip, cloud_cover, wind, fog, Site, .direction = "down") %>%
  ungroup() %>%
  # add 1s and 2s to indicate whether an individual was resighted 
  mutate(sighted = ifelse(is.na(band_as_seen), 1, 2)) %>%
  # merge with the nesting peruiod data
  left_join(
    nests.br.periods,
    by = join_by("band_id" == "male_id", "year" == "year"),
    relationship = "many-to-many") %>%
  # Check if each female is in any of the breeding states 
  mutate(in_state = !is.na(start_time) &
           !is.na(end_time) &
           date >= start_time &
           date < end_time) %>%
  # If a bird is in none of the states, we will give it the "preparation" state
  group_by(band_id, date) %>%
  filter(in_state == T | all(in_state == F)) %>%
  mutate(state = ifelse(in_state == F, "preparation", state)) %>%
  distinct(band_id, date, .keep_all = T) %>%
  # then we can define pre-breeding and post-breeding states 
  group_by(band_id, year) %>%
  mutate(breeding.start = min(start_time),
         breeding.end = max(end_time),
         state = case_when(!is.na(breeding.start) & date <= breeding.start ~ "prebreeding",
                           !is.na(breeding.end) & date > breeding.end ~ "postbreeding",
                           is.na(breeding.start) & yday(date) <= 193 ~ "prebreeding",
                           is.na(breeding.end) & yday(date) > 227 ~ "postbreeding",
                           is.na(breeding.end) & yday(date) > 227 & yday(date) < 227 ~ "preparation",
                           .default = state)) %>%
  arrange(band_id, date) %>%
  mutate(period = factor(period, levels = c("prebreeding", "breeding", "postbreeding")),
         state = factor(state, levels = c("prebreeding", "preparation", "nestling care", "fledgling care", "postbreeding"))) %>%
  # lastly, we will modify the "start-time" and "end_time" column so they include updated information 
  group_by(band_id) %>%
  mutate(state_id = cumsum(state != lag(state, default = first(state)))) %>%
  group_by(band_id, state_id) %>%
  mutate(start_time = min(date),
         end_time = max(date)) %>%
  ungroup()

# Merge the male and female resighting datasets with states 
br.survey.states <- bind_rows(br.survey.states.m, br.survey.states.f)

# Plot resightings with the old periods ----
male_res_per <- ggplot() + 
  geom_tile(data = br.survey.states.m, aes(x= as.factor(date), y = likely_band, fill = period), colour = "white") +
  scale_fill_manual(values = c("breeding" = "#E69F00", "prebreeding" ="#56B4E9", "postbreeding" = "#009E73"), 
                           labels = c("breeding" = "Breeding", "prebreeding" = "Prebreeding",
                                      "postbreeding" = "Postbreeding"))+
  labs(fill = "Period",
      title  = "(a) Male resightings") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 11))

female_res_per <- ggplot() + 
  geom_tile(data = br.survey.states.f, aes(x= as.factor(date), y = likely_band, fill = period), colour = "white") + 
  scale_fill_manual(values = c("breeding" = "#E69F00", "prebreeding" ="#56B4E9", "postbreeding" = "#009E73"), 
                    labels = c("breeding" = "Breeding", "prebreeding" = "Prebreeding",
                               "postbreeding" = "Postbreeding"))+
  
  labs(fill = "Period",
       title  = "(b) Female resightings") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 11))

res_periods <- male_res_per / female_res_per

ggsave(plot = res_periods , device = "png" ,filename = "Resightings_period.png", path = fig.dir, dpi=300, units = "cm", width = 25, height = 10)

# Save the the processed data ---
write_csv(br.survey.states, paste0(local.data,"/breeding_survival/SAVS_resightings_pro.csv"))

# Plot resightings with new states ----
male_res_sta <- ggplot() + 
  geom_tile(data = br.survey.states.m, aes(x= as.factor(date), y = likely_band, fill = state), colour = "white") + 
  scale_fill_discrete(labels = c("breeding" = "Breeding", "prebreeding" = "Prebreeding",
                               "postbreeding" = "Postbreeding", "fledgling care" = "Fledgling care",
                               "nestling care" = "Nestling care",
                               "preparation" = "Incubation\n & nest bulding"))+
  labs(fill = "Period",
       title  = "(a) Male resightings") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 11))

female_res_sta <- ggplot() + 
  geom_tile(data = br.survey.states.f, aes(x= as.factor(date), y = likely_band, fill = state), colour = "white") + 
  scale_fill_discrete(labels = c("breeding" = "Breeding", "prebreeding" = "Prebreeding",
                                 "postbreeding" = "Postbreeding", "fledgling care" = "Fledgling care",
                                 "nestling care" = "Nestling care",
                                 "preparation" = "Incubation\n & nest bulding"))+
  labs(fill = "Period",
       title  = "(b) Female resightings") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 11))

res_states<- male_res_sta / female_res_sta

ggsave(plot = res_states , device = "png" ,filename = "Resightings_states.png", path = fig.dir, dpi=300, units = "cm", width = 25, height = 10)

# Save the the processed data ---
write_csv(br.survey.states, paste0(local.data,"/breeding_survival/SAVS_resightings_pro.csv"))



