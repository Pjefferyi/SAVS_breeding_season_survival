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
detec_hist  <- read_csv(paste0(local.data,"/breeding_survival/SAVS_resightings_pro.csv"))

## Complete the data by adding non-detections as NA ----
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
delta.days <- detec_hist %>% ungroup() %>% distinct(date) %>%
  arrange(date) %>% 
  mutate(time.diff = as.numeric(difftime(date, lag(date), units = "days"))) %>%
  filter(!is.na(time.diff )) %>%
  pull(time.diff)

# For modelling, obtain a list of sexes ----
sex <- detec_hist %>% distinct(likely_band, sex) %>% filter(!is.na(sex)) %>%
  mutate(sex = ifelse(sex == "M", 1, 2)) %>% pull(sex)

# Survival model accounting for time interval between sampling efforts ----

# dt[t] represents the length of time between occasion t and t+1
# phi_per_unit is your baseline survival rate per unit time (e.g., per year)

code <- nimbleCode({
  # Priors
  phi_per_unit ~ dunif(0, 1)
  for (t in 1:(T-1)) {
    p[t] ~ dunif(0, 1)          # Capture probability at each occasion
  }
  
  # Process & Observation models
  for (i in 1:N) {
    # Define state at first capture
    z[i, first[i]] <- 1 
    
    for (t in (first[i]):(T-1)) {
      # Calculate survival for this specific interval length
      phi_interval[i, t] <- phi_per_unit^dt[t] 
      
      # State process: Survival
      z[i, t+1] ~ dbern(z[i, t] * phi_interval[i, t])
      
      # Observation process: Resighting
      y[i, t+1] ~ dbern(z[i, t+1] * p[t])
      
    }
  }
})

# We get the occasion of first capture for all individuals 
first <- apply(y, 1, function(x) min(which(x ==2)))

# Now we specify the constants 
my.constants <- list(N = nrow(y),   # number of animals
                     T = ncol(y),   # number of sampling occasions
                     first = first) # first resighting occasion for each individual 

# Put data in a list
my.data <- list(y = y)

# Now we write a function for the initial values 
# here, for latent states, we assume that individuals are alive through the study period 
zinits <- y # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() {
  list(
    phi = matrix(
      runif(2 * (my.constants$T - 1), 0.5, 0.99),
      nrow = 2,
      ncol = my.constants$T - 1
    ),
    p = matrix(
      runif(2 * (my.constants$T - 1), 0, 1),
      nrow = 2,
      ncol = my.constants$T - 1
    ),
    z = zinits
  )
}

# We specify the parameter sto monitor
parameters.to.save <- c("phi", "p")

# set MCM details 
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2




hmm.phisexpt.ni <- nimbleCode({
  
  # Initial state probabilities
  delta[1] <- 1
  delta[2] <- 0
  
  # Interval-specific survival by sex
  for (s in 1:2) {
    for (t in 1:(T-1)) {
      
      phi[s,t] ~ dunif(0, 1)
      
      gamma[1,1,t,s] <- phi[s,t]
      gamma[1,2,t,s] <- 1 - phi[s,t]
      gamma[2,1,t,s] <- 0
      gamma[2,2,t,s] <- 1
      
      # Detection probability
      p[s,t] ~ dunif(0, 1)
      
      omega[1,1,t,s] <- 1 - p[s,t]
      omega[1,2,t,s] <- p[s,t]
      omega[2,1,t,s] <- 1
      omega[2,2,t,s] <- 0
    }
  }
  
  # Likelihood
  for (i in 1:N) {
    
    z[i,first[i]] ~ dcat(delta[1:2])
    
    for (j in (first[i]+1):T) {
      
      # State transition using the bird's sex
      z[i,j] ~ dcat(
        gamma[z[i,j-1], 1:2, j-1, sex[i]]
      )
      
      # Observation model
      y[i,j] ~ dcat(
        omega[z[i,j], 1:2, j-1, sex[i]]
      )
    }
  }
})

# We get the occasion of first capture for all individuals 
first <- apply(y, 1, function(x) min(which(x ==2)))

# Now we specify the constants 
my.constants <- list(N = nrow(y),   # number of animals
                     T = ncol(y),   # number of sampling occasions
                     first = first, # first resighting occasion for each individual 
                     sex = sex)

# Put data in a list
my.data <- list(y = y)

# Now we write a function for the initial values 
# here, for latent states, we assume that individuals are alive through the study period 
zinits <- y # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() {
  list(
    phi = matrix(
      runif(2 * (my.constants$T - 1), 0.5, 0.99),
      nrow = 2,
      ncol = my.constants$T - 1
    ),
    p = matrix(
      runif(2 * (my.constants$T - 1), 0, 1),
      nrow = 2,
      ncol = my.constants$T - 1
    ),
    z = zinits
  )
}

# We specify the parameter sto monitor
parameters.to.save <- c("phi", "p")

# set MCM details 
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2