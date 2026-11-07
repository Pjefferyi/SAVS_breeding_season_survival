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

# We will start by building a simpleCormack Jolly Seber model for Kent Island Savannah sparrows
# We will not have distinct states for the breeding status of individuals 
# But we will have time varying survival and detection probabilities 
# For now we will assume that the interval between resighting efforts is constant, even though it is not always 

# Set local repository in which to store data ----
local_data <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/PhD_Thesis/Data"

# Load processed SAVS resighting data 
detec_hist  <- read_csv(paste0(local.data,"/breeding_survival/SAVS_resightings_pro.csv"))

## Complete the data by adding non-detections as NA ----
detec_hist <- detec_hist %>% complete(date, likely_band, fill = list(period = "Not seen")) %>%
  group_by(date)%>%
  fill(surveyor, start_time, end_time, precip, cloud_cover, wind, fog, Site, .direction = "downup") %>%
  ungroup() %>%
  # add 1s and 2s to indicate whether an individual was resighted 
  mutate(sighted = ifelse(is.na(band_as_seen), 1, 2))

# We convert our data to a resighting matrix ---- 
y <- detec_hist %>% select(date, sighted, likely_band) %>% group_by(date) %>% 
  distinct(likely_band, .keep_all = T) %>%
  pivot_wider(names_from = date, values_from = sighted, id_cols = likely_band)

y <- y %>% select(-likely_band)%>% as.matrix()

# Survival Model with time-varying survival and recapture probability 
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
                     first = first) # occasion of first capture for each individual

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
         params = "p")

# You can visualize the traceplot
MCMCtrace(object = mcmc.phitpt,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "phi")

#  We can create a dataframe of the results and plot survival
df.surv.prob <- data_frame("days" = unique(detec_hist$date)[2:61], "survival.prob" = mcmc.sm$mean[1:60])

ggplot(data = df.surv.prob, aes(x = days, y = survival.prob)) + 
  geom_path() + 
  theme_bw()

# Model with constant parameters ----
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


# Re-design the tiem varying model to account for the length of time intervals ----

# We have to get the time intervals between resightign events
delta.days <- detec_hist %>% ungroup() %>% distinct(date) %>%
  arrange(date) %>% 
  mutate(time.diff = as.numeric(difftime(date, lag(date), units = "days"))) %>%
  filter(!is.na(time.diff )) %>%
  pull(time.diff)

# Survival Model with time-varying survival and recapture probability 
hmm.phitpt <- nimbleCode({
  
  # initial parameter values 
  delta[1] <- 1                 # Pr(alive t = 1) = 1
  delta[2] <- 0                 # Pr(dead t = 1) = 0
  
  # Daily survival probability
  S ~ dunif(0, 1)
  
  for (t in 1:(T-1)){
    phi[t] <- pow(S, delta.days[t]) # Survival over the entire time interval 
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1           # Pr(dead t -> dead t+1)
    p[t] ~ dunif(0, 1)          # prior detection at occasion t
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
                     first = first, # first resighting occasion for each individual 
                     delta.days = delta.days)

# Put data in a list
my.data <- list(y = y)

# Now we write a function for the initial values 
# here, for latent states, we assume that individuals are alive through the study period 
zinits <- y # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(S = runif(1, 0.95, 1),
                                  #phi = runif(my.constants$T-1,0,1),
                                  p = runif(my.constants$T-1,0,1),
                                  z = zinits)

# We specify the parameter sto monitor
parameters.to.save <- c("S", "phi", "p")

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

#  We can create a dataframe of the results and plot survival
df.surv.prob <- data_frame("days" = unique(detec_hist$date)[1:60], "survival.prob" = mcmc.sm$mean[62:121])

ggplot(data = df.surv.prob, aes(x = days, y = survival.prob)) + 
  geom_path() + 
  theme_bw()

