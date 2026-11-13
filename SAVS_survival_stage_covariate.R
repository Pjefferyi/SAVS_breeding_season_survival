# SAVS breeding survival model 

# Resources 

# Olivier Gimenez guide: https://oliviergimenez.github.io/bayesian-hmm-cr-workshop-valencia/worksheets/2_demo.html 
# Olivier Giminez guide2: https://oliviergimenez.github.io/banana-book/intronimble.html 
# Introduction to Nimble: https://cran.r-project.org/web/packages/nimbleEcology/vignettes/Introduction_to_nimbleEcology.html 

# Load libraries 
library(nimbleEcology)
library(nimble)
library(tidyverse)
library(MCMCvis)
library(anytime)

# Set a seed for reproducibility 
mySeed <- 0
set.seed(mySeed)

# We will start by building a simpleCormack Jolly Seber model for Kent Island Savannah sparrows
# We will not have distinct states for the breeding status of individuals 
# But we will have time varying survival and detection probabilities 
# For now we will assume that the interval between resighting efforts is constant, even though it is not always 

# Set local repository in which to store data ----
local_data <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/PhD_Thesis/Data"

# Load processed SAVS resighting data 
detec_hist_full <- read_csv(paste0(local_data,"/breeding_survival/SAVS_resightings_pro.csv"))

# We can break this history down in a full (all days) and survey only dataset (only survey days) 
detec_hist <- detec_hist_full %>% filter(!is.na(period))

## Complete the survey data by adding non-detections as NA ----
detec_hist_fill <- detec_hist %>% complete(date, nesting(likely_band), fill = list(period = "Not seen")) %>%
  group_by(date)%>%
  fill(surveyor, start_time, end_time, precip, cloud_cover, wind, fog, Site, .direction = "down") %>%
  ungroup() %>%
  # add 1s and 2s to indicate whether an individual was resighted 
  mutate(sighted = ifelse(is.na(band_as_seen), 1, 2))

# We convert our data to a resighting matrix ---- 
y <- detec_hist_fill %>% select(date, sighted, likely_band) %>% group_by(date) %>% 
  distinct(likely_band, .keep_all = T) %>%
  pivot_wider(names_from = date, values_from = sighted, id_cols = c(likely_band)) %>%
  mutate(likely_band = factor(likely_band, levels = unique(detec_hist$likely_band))) %>%
  arrange(likely_band)

y <- y %>% select(-likely_band)%>% as.matrix()

# For modelling, obtain the time between resighting attempts ----
interval.info <-detec_hist %>% ungroup() %>% distinct(date) %>%
  arrange(date) %>% 
  mutate(time.diff = as.numeric(difftime(date, lag(date), units = "days")),
         interval.type = case_when(time.diff < 10 ~ 1, 
                                   time.diff > 10 & time.diff < 60 ~ 2,
                                   time.diff > 60 ~ 3))
delta.days <- interval.info %>%
  filter(!is.na(time.diff )) %>%
  pull(time.diff)

# For modelling, define intervals as April - August (breeding state timing known), August - September (post breeding period), or September-April (mix of migration and post-breeding, timings unknown)
interval.type <- interval.info %>% filter(!is.na(time.diff )) %>%
  pull(interval.type)

# For modelling, define the amount of time spent in each state ----
days.state <- array(
  0,
  dim = c(122, nocc - 1, 4)
)


# For modelling, obtain a list of sexes ----
sex <- detec_hist %>% distinct(likely_band, sex) %>% filter(!is.na(sex)) %>%
  mutate(sex = ifelse(sex == "M", 1, 2)) %>% pull(sex)
