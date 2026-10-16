# Generate statistics for breeding success on Kent Island

## load packages ---
library(tidyverse)
library(ggplot2)
library(anytime)

Sys.setlocale("LC_ALL", "C")

## Path to Jelany's local data ---
local.data <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/PhD_Thesis/Data"

## Load Savannah sparrow nest data 

## load SAVS dataset 
individs <- read_csv(paste0(local.data,"/Kent_morpho_nest_ind/SAVS_dataset_backup_2026_08_03/JD_individuals_20260803.csv"),
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
                     )) %>% mutate(banding_date= anytime(ymd(bandign_date)))

morph <- read_csv(paste0(local.data,"/Kent_morpho_nest_ind/SAVS_dataset_backup_2026_08_03/JD_morphology_20260803.csv"),
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

nests <- read_csv(paste0(local.data,"/Kent_morpho_nest_ind/SAVS_dataset_backup_2026_08_03/JD_nests_20260803.csv"),
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
                  )) %>% mutate(found_date = anytime(ymd(found_date)))


# Calculate plot the frequency of nests in each category for 2026 and 2025
nests_rec <- nests %>% filter(year %in% c(2023, 2024, 2025, 2026))

ggplot() +
  geom_bar(data = nests_rec, aes(fate, group = factor(year), fill = factor(year)),
           position = "dodge")

nests_rec %>% filter(year == 2026) %>% group_by(fate) %>% tally()

# plot nest fail date over time in  2026
nests_2026 <- nests %>% filter(year %in% c(2026))

ggplot() +
  geom_bar(data = nests_2026, aes(found_date, group = factor(fate), fill = factor(fate)),
           position = "dodge")


