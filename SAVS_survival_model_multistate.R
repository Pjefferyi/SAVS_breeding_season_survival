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

## Complete the data by adding non-detections as NA ----
detec_hist_fill <- detec_hist %>% complete(date, nesting(likely_band), fill = list(period = "Not seen")) %>%
  group_by(date)%>%
  fill(surveyor, start_time, end_time, precip, cloud_cover, wind, fog, Site, .direction = "down") %>%
  ungroup() %>%
  # add 1s and 2s to indicate whether an individual was resighted 
  mutate(sighted = case_when(is.na(band_as_seen) ~ 0,
                             !is.na(band_as_seen) & state == "prebreeding" ~1,
                             !is.na(band_as_seen) & state == "preparation" ~2,
                             !is.na(band_as_seen) & state == "nestling care" ~3,
                             !is.na(band_as_seen) & state =="fledgling care" ~4,
                             !is.na(band_as_seen) & state == "postbreeding" ~5,))

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

# Multi-state CJS model ----

multistate <- nimbleCode({
  
  # Priors
  
  # survival 
  
  phi1 ~ dunif(0, 1) # pre-breeding survival
  phi2 ~ dunif(0, 1) # preparation survival
  phi3 ~ dunif(0, 1) # feeding nestlings survival
  phi4 ~ dunif(0, 1) # feeding fledglings survival
  phi5 ~ dunif(0, 1) # post-breeding survival 
  
  # # Regular transition probabilities
  # 
  # psi12 ~ ddirch(alpha[1:5]) # transition from prebreeding season to preparation  
  # psi23 ~ ddirch(alpha[1:5]) # transition from preparation season to feeding nestlings
  # psi34 ~ ddirch(alpha[1:5]) # transition from feeding nestlings season to feeding fledglings
  # psi45 ~ ddirch(alpha[1:5]) # transition from feeding fledglings season to post-breeding 
  # psi51 ~ ddirch(alpha[1:5]) # transition from post-breeding to prebreeding 
  # 
  # # Irregular transition probabilities (generally for nest failures, floater) 
  # 
  # psi32 ~ ddirch(alpha[1:5]) # transition from feeding nestlings to preparation  
  # psi42 ~ ddirch(alpha[1:5]) # transition from feeding fledglings to preparation 
  # psi35 ~ ddirch(alpha[1:5]) # transition from feeding nestlings to postbreeding
  # psi25 ~ ddirch(alpha[1:5]) # transition from preparation to postbreeding 
  
  # Transition probabilities conditional on survival
  
  # Pre-breeding: remain or move to preparation
  psi1[1:2] ~ ddirch(alpha[1:2])
  
  # Preparation: remain, progress, or move to post-breeding
  psi2[1:3] ~ ddirch(alpha[1:3])
  
  # Nestling care: move to preparation, remain, progress, or post-breeding
  psi3[1:4] ~ ddirch(alpha[1:4])
  
  # Fledgling care: move to preparation, remain, or post-breeding
  psi4[1:3] ~ ddirch(alpha[1:3])
  
  # Post-breeding: remain or move to pre-breeding
  psi5[1:2] ~ ddirch(alpha[1:2])
  
  # Resighting probability (numbers here follow survival probabilities)
  
  p1 ~ dunif(0, 1)
  p2 ~ dunif(0, 1)
  p3 ~ dunif(0, 1)
  p4 ~ dunif(0, 1)
  p5 ~ dunif(0, 1)
  
  # Initial state probabilities
  # 
  # delta[1] <- pi1          # Pr(alive in 1 at t = first)
  # delta[2] <- 1 - pi1      # Pr(alive in 2 at t = first)
  # delta[3] <- 0            # Pr(dead at t = first) = 0
  
  # State transition matrix (grouped)
  # probabilities of state z(t+1) given z(t)
  
  # #prebreeding
  # gamma[1,1] <- phi1 * (1 - psi12)
  # gamma[1,2] <- phi1 * psi12
  # gamma[1,3] <- 0
  # gamma[1,4] <- 0
  # gamma[1,5] <- 0
  # gamma[1,6] <- 1- phi1
  # 
  # # preparation
  # gamma[2,1] <- 0
  # gamma[2,2] <- phi2 * (1 - psi23 - psi25)
  # gamma[2,3] <- phi2 * psi23
  # gamma[2,4] <- 0
  # gamma[2,5] <- phi2 * psi25
  # gamma[2,6] <- 1 - phi2
  # 
  # # nestling care
  # gamma[3,1] <- 0
  # gamma[3,2] <- phi3 * psi32
  # gamma[3,3] <- phi3 * (1 - psi32 - psi34 - psi35)
  # gamma[3,4] <- phi3 * psi34
  # gamma[3,5] <- phi3 * psi35
  # gamma[3,6] <- 1 - phi3
  # 
  # # Fledgeling care
  # gamma[4,1] <- 0
  # gamma[4,2] <- phi4 * psi42
  # gamma[4,3] <- 0
  # gamma[4,4] <- phi4 * (1 - psi42 - psi45)
  # gamma[4,5] <- phi4 * psi45
  # gamma[4,6] <- 1-phi4
  # 
  # # postbreeding 
  # gamma[5,1] <- phi5 * psi51
  # gamma[5,2] <- 0
  # gamma[5,3] <- 0
  # gamma[5,4] <- 0
  # gamma[5,5] <- phi5 * (1 - psi51)
  # gamma[5,6] <- 1 - phi5
  # 
  # # Dead state
  # gamma[6,1] <- 0
  # gamma[6,2] <- 0
  # gamma[6,3] <- 0
  # gamma[6,4] <- 0
  # gamma[6,5] <- 0
  # gamma[6,6] <- 1
  
  # Pre-breeding
  gamma[1,1] <- phi1 * psi1[1]
  gamma[1,2] <- phi1 * psi1[2]
  gamma[1,3] <- 0
  gamma[1,4] <- 0
  gamma[1,5] <- 0
  gamma[1,6] <- 1 - phi1
  
  # Preparation
  gamma[2,1] <- 0
  gamma[2,2] <- phi2 * psi2[1]
  gamma[2,3] <- phi2 * psi2[2]
  gamma[2,4] <- 0
  gamma[2,5] <- phi2 * psi2[3]
  gamma[2,6] <- 1 - phi2
  
  # Nestling care
  gamma[3,1] <- 0
  gamma[3,2] <- phi3 * psi3[1]
  gamma[3,3] <- phi3 * psi3[2]
  gamma[3,4] <- phi3 * psi3[3]
  gamma[3,5] <- phi3 * psi3[4]
  gamma[3,6] <- 1 - phi3
  
  # Fledgling care
  gamma[4,1] <- 0
  gamma[4,2] <- phi4 * psi4[1]
  gamma[4,3] <- 0
  gamma[4,4] <- phi4 * psi4[2]
  gamma[4,5] <- phi4 * psi4[3]
  gamma[4,6] <- 1 - phi4
  
  # Post-breeding
  gamma[5,1] <- phi5 * psi5[2]
  gamma[5,2] <- 0
  gamma[5,3] <- 0
  gamma[5,4] <- 0
  gamma[5,5] <- phi5 * psi5[1]
  gamma[5,6] <- 1 - phi5
  
  # Dead
  gamma[6,1] <- 0
  gamma[6,2] <- 0
  gamma[6,3] <- 0
  gamma[6,4] <- 0
  gamma[6,5] <- 0
  gamma[6,6] <- 1
  
  # probabilities of y(t) given z(t)
  # (read as omega[y(t),z(t)] = omega[Observation,State])
  
  omega[1,1] <- 1 - p1     
  omega[1,2] <- p1         
  omega[1,3] <- 0  
  omega[1,4] <- 0 
  omega[1,5] <- 0 
  omega[1,6] <- 0  
  
  omega[2,1] <- 1 - p2     
  omega[2,2] <- 0         
  omega[2,3] <- p2  
  omega[2,4] <- 0  
  omega[2,5] <- 0  
  omega[2,6] <- 0  
  
  omega[3,1] <- 1 - p3     
  omega[3,2] <- 0         
  omega[3,3] <- 0  
  omega[3,4] <- p3  
  omega[3,5] <- 0  
  omega[3,6] <- 0    
  
  omega[4,1] <- 1 - p4     
  omega[4,2] <- 0         
  omega[4,3] <- 0  
  omega[4,4] <- 0  
  omega[4,5] <- p4  
  omega[4,6] <- 0   
  
  omega[5,1] <- 1 - p5     
  omega[5,2] <- 0         
  omega[5,3] <- 0  
  omega[5,4] <- 0  
  omega[5,5] <- 0 
  omega[5,6] <- p5   
  
  omega[6,1] <- 1     
  omega[6,2] <- 0         
  omega[6,3] <- 0  
  omega[6,4] <- 0  
  omega[6,5] <- 0  
  omega[6,6] <- 0   
  
  # likelihood
  for (i in 1:N){
    # latent state at first capture
    z[i,first[i]] <- y[i,first[i]] - 1
    for (t in (first[i]+1):K){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1],1:6])
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:6])
    }
  }
})

# occasions of first capture
first <- apply(y, 1, function(x) min(which(x !=0)))

# constants
my.constants <- list(first = first, 
                     K = ncol(y), 
                     N = nrow(y),
                     alpha = c(1, 1, 1, 1))
# data
my.data <- list(y = y + 1)

# Generate Dirichlet random values
rdirichlet <- function(n, alpha) {
  x <- matrix(
    rgamma(n * length(alpha),
           shape = alpha,
           rate = 1),
    nrow = n
  )
  x / rowSums(x)
}


# initial values 
zinits <- y

# non-detections become alive in site A or B
zinits[zinits==0] <- sample(c(1,2,3,4,5), sum(zinits==0), replace = TRUE) 
initial.values <- function(){list(phi1 = runif(1, 0, 1), 
                                  phi2 = runif(1, 0, 1),
                                  phi3 = runif(1, 0, 1), 
                                  phi4 = runif(1, 0, 1),
                                  phi5 = runif(1, 0, 1), 
                                  
                                  p1 = runif(1, 0, 1),
                                  p2 = runif(1, 0, 1), 
                                  p3 = runif(1, 0, 1),
                                  p4 = runif(1, 0, 1), 
                                  p5 = runif(1, 0, 1),
                                  
                                  # psi12 = runif(1, 0, 1), 
                                  # psi23 = runif(1, 0, 1), 
                                  # psi34 = runif(1, 0, 1), 
                                  # psi45 = runif(1, 0, 1), 
                                  # psi51 = runif(1, 0, 1), 
                                  # psi32 = runif(1, 0, 1), 
                                  # psi42 = runif(1, 0, 1), 
                                  # psi35 = runif(1, 0, 1), 
                                  # psi25 = runif(1, 0, 1), 
                                  
                                  psi1 = as.numeric(rdirichlet(1, alpha[1:2])),
                                  psi2 = as.numeric(rdirichlet(1, alpha[1:3])),
                                  psi3 = as.numeric(rdirichlet(1, alpha[1:4])),
                                  psi4 = as.numeric(rdirichlet(1, alpha[1:3])),
                                  psi5 = as.numeric(rdirichlet(1, alpha[1:2])),
                                  
                              
                                  z = zinits)}
# parameters to monitor
parameters.to.save <- c("phi1", "phi2",
                        "phi3", "phi4",
                        "phi5",
                        
                        # "psi12", "psi23",
                        # "psi34", "psi45",
                        # "psi51", "psi32",
                        # "psi42", "psi35",
                        # "psi25",
                        
                        "psi1", "psi2",
                        "psi3", "psi4",
                        "psi5",
                         
                        "p1", "p2",
                        "p3", "p4",
                        "p5")

# MCMC details
n.iter <- 40000
n.burnin <- 1000
n.chains <- 2

mcmc.multi.state = nimbleMCMC(code = multistate, 
                              constants = my.constants,
                              data = my.data,              
                              inits = initial.values,
                              monitors = parameters.to.save,
                              niter = n.iter,
                              nburnin = n.burnin, 
                              nchains = n.chains,
                              dimensions = list(gamma = c(6, 6)))
