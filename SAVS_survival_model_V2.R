# SAVS breeding survival model 

# Resources 

# Olivier Gimenez guide: https://oliviergimenez.github.io/bayesian-hmm-cr-workshop-valencia/worksheets/2_demo.html 
# Olivier Giminez guide2 : https://oliviergimenez.github.io/banana-book/intronimble.html 
# Introduction to Nimble: https://cran.r-project.org/web/packages/nimbleEcology/vignettes/Introduction_to_nimbleEcology.html 

# Load libraries 
library(nimbleEcology)
library(nimble)
library(tidyverse)
library(MCMCvis)
library(anytime)


# We will start by building a simpleCormack Jolly Seber model for Kent Island Savannah sparrows
# We will not have distinct states for the breeding status of individuals 
# But we will have time varying survival and detection probabilities 
# For now we will assume that the interval between resighting efforts is constant, even though it is not always 

# Set local repository in which to store data ----
local_data <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/PhD_Thesis/Data"

# First, let's load our Savannah Sparrow data ----

## load individual SAVS data ----
individs <- read_csv(paste0(local_data,"/Kent_morpho_nest_ind/SAVS_data_backup_20251113/JD_individuals_20251113.csv"),
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
                       banding_date = col_character(),
                       natal_nestkey = col_character(),
                       mother_color = col_character(),
                       father_color = col_character(),
                       Notes = col_character()
                     )) %>% mutate(banding_date = anytime(ymd(banding_date)))

## Load the resighting data ----
br.survey <- read_csv(paste0(local_data,"/breeding_survival/SAVS_resightings_20250930.csv")) %>%
  mutate(date = anytime(dmy(date))) %>% filter(!(likely_band %in% c("unknown", ".A", "A", "W", ".W", "W."))) %>%
  mutate(period = case_when(period == "breeding" ~ "Breeding",
                            period == "prebreeding" ~ "Pre-breeding"))

## Load morphology data to get data from the fall banding period ----
morph <- read_csv(paste0(local_data,"/Kent_morpho_nest_ind/SAVS_data_backup_20251029/JD_morphology_20251029.csv"),
                  col_types = cols(
                    band_id = col_character(),
                    location = col_character(),
                    weight = col_double(),
                    tarsus = col_double(),
                    wing_l = col_double(),
                    tail_l = col_double(),
                    bill_l = col_double(),
                    bill_d = col_double(),
                    date_measured = col_character(),
                    time_measured = col_double(),
                    sex = col_character(),
                    age_in_yrs = col_double(),
                    age_in_d = col_double(),
                    color = col_character(),
                    year = col_double(),
                    fat_score = col_double(),
                    bander_id = col_character(),
                    Notes = col_character()
                  ))

## Take morphology data from the fall banding period and merge it with the resighting data ----
fall_bd <- morph %>% mutate(date_measured = anytime(ymd(date_measured))) %>%
  filter(date_measured > "2025-08-29" & date_measured < "2025-09-08") %>%
  group_by(location) %>%
  distinct(color, .keep_all = T) %>% filter(color %in% br.survey$likely_band) %>% 
  group_by(location, date_measured) %>%
  summarize(likely_band = unique(color)) %>% mutate(period = "Post-breeding") %>%
  rename(Location_first_seen = location,
         date = date_measured) %>%
  mutate(surveyor = "Fall banding",
         start_time = NA,
         end_time = NA,
         precip = NA,
         cloud_cover = NA, 
         fog = NA, 
         wind = NA,
         band_as_seen = NA,
         Site = "S")

## Complete the data by adding non-detections as NA ----
detec_hist <- rbind(br.survey, fall_bd)
detec_hist <- detec_hist %>% complete(date, likely_band, fill = list(period = "Not seen")) %>%
  filter(!(likely_band %in% c("A.", "YR.B", "SL.R", "WW.R", "WN.N"))) %>%
  group_by(date)%>%
  fill(surveyor, start_time, end_time, precip, cloud_cover, wind, fog, Site, .direction = "downup") %>%
  ungroup() %>%
  # add 1s and 2s to indicate whether an individual was resighted 
  mutate(sighted = ifelse(is.na(band_as_seen), 1, 2))

# Now add individual characteristics from the individual dataset ----
individs.prep <- individs %>% filter(adult_color %in% detec_hist$likely_band) %>%
  group_by(adult_color) %>% filter(banding_date == max(banding_date))
detec_hist <- merge(x = detec_hist, y = individs.prep, by.x = "likely_band", by.y = "adult_color", all = T)
  # We also convert resighting occasions to 1 and 0s 

# Calculate the age of individuals ----
detec_hist <- detec_hist %>% mutate(age = 2025 - year(banding_date),
                                    age_category = ifelse(age > 1, "ASY","SY")) %>%
  arrange(band_id, date)

# We convert our data to a resighting matrix ---- 
y <- detec_hist %>% select(date, sighted, band_id) %>% group_by(date) %>% 
  distinct(band_id, .keep_all = T) %>%
  pivot_wider(names_from = date, values_from = sighted, id_cols = band_id)

y <- y %>% select(-band_id)%>% as.matrix()
  
# Now we have to write a model 
hmm.phitpt <- nimbleCode({
  # parameters
  delta[1] <- 1                 # Pr(alive t = 1) = 1
  delta[2] <- 0                 # Pr(dead t = 1) = 0
  for (t in 1:(T-1)){
    phi[t] ~ dunif(0, 1)        # prior survival
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1           # Pr(dead t -> dead t+1)
    p[t] ~ dunif(0, 1)          # prior detection
    omega[1,1,t] <- 1 - p[t]    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t]        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0           # Pr(dead t -> detected t)
  }
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
    }
  }
})

# We get the occasion of first capture for all individuals 
first <- apply(y, 1, function(x) min(which(x ==2)))

# Now we specify the constants 
my.constants <- list(N = nrow(y),   # number of animals
                     T = ncol(y),   # number of sampling occasions
                     first = first)

# Put data in a list
my.data <- list(y = y)

# Now we write a function for the initial values 
# here, for latent states, we assume that individuals are alive through the study period 
zinits <- y # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(my.constants$T-1,0,1),
                                  p = runif(my.constants$T-1,0,1),
                                  z = zinits)

# We specify the parameter sto monitor
parameters.to.save <- c("phi", "p")

# set MCM details 
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

# And we can run Nimble 
mcmc.phitpt <- nimbleMCMC(code = hmm.phitpt,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)

# summarize and visualize the MCMC output with the package MCMCvis
mcmc.sm <- MCMCsummary(object = mcmc.phitpt, round = 2)
mcmc.sm$param <- row.names(mcmc.sm)

# We can visualize the posterior distribution of thetha 
MCMCplot(object = mcmc.phitpt, 
         params = "phi")

# You can visualize the traceplot
MCMCtrace(object = mcmc.phitpt,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "phi")

# Now let's try a model with constant parameters
hmm.phip <- nimbleCode({
  phi ~ dunif(0, 1)      # prior survival
  p ~ dunif(0, 1)        # prior detection
  # likelihood
  gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0        # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1        # Pr(dead t -> dead t+1)
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  omega[1,1] <- 1 - p    # Pr(alive t -> non-detected t)
  omega[1,2] <- p        # Pr(alive t -> detected t)
  omega[2,1] <- 1        # Pr(dead t -> non-detected t)
  omega[2,2] <- 0        # Pr(dead t -> detected t)
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# And we can run Nimble 
hmm.phip <- nimbleMCMC(code = hmm.phip,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)
