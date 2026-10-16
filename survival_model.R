# Estimating survival rates during the breeding period on Kent island

# Resources 

# Olivier Gimenez guide: https://oliviergimenez.github.io/bayesian-hmm-cr-workshop-valencia/worksheets/2_demo.html 
# Olivier Giminez guide2 : https://oliviergimenez.github.io/banana-book/intronimble.html 
# Introduction to Nimble: https://cran.r-project.org/web/packages/nimbleEcology/vignettes/Introduction_to_nimbleEcology.html 

# Load libraries 
library(nimbleEcology)
library(nimble)
library(tidyverse)
library(MCMCvis)

# We start by building a very simple surival model 
model <-  nimbleCode({
  # likelihood
  survived ~ dbinom(theta, released)
  # prior
  theta ~ dunif(0, 1)
  # derived quantity
  lifespan <- -1/log(theta)
})

# read in data 
my.data <- list(released = 57, survived = 19)

# identify parameters that NIMBLE should keep track of 
parameters.to.save <- c("theta", "lifespan")

# Specify intial values for 3 chains (We could also have generated some random variables, remember to set a seed)
init1 <- list(theta = 0.1)
init2 <- list(theta = 0.5)
init3 <- list(theta = 0.9)
initial.values <- list(init1, init2, init3)

# determine how many chains we would like to run, on two many chains
n.iter <- 5000
n.burnin <- 1000
n.chains <- 3

# We can now run our model 
# You probably want to set a seed using setSeed here for reproducibility 
mcmc.output <- nimbleMCMC(code = model,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)

# We can visualize the results of the theta argument (i.e. the probability of survival)
mcmc.output$chain1[,"theta"] %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = value), color = "white") +
  labs(x = "survival probability")

# You can also summarize and visualize the MCMC output with the package MCMCvis
MCMCsummary(object = mcmc.output, round = 2)

# We can visualize the posterior distribution of thetha 
MCMCplot(object = mcmc.output, 
         params = 'theta')

# You can visualize the traceplot
MCMCtrace(object = mcmc.output,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "theta")

# You can add some divergence diagnostics
MCMCtrace(object = mcmc.output,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE, # add Rhat
          n.eff = TRUE, # add eff sample size
          params = "theta")

# Now we run the Dipper example 

#Load Dipper data 
dipper <- read_csv("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/PhD_Thesis/Chapter_3_breeding_survival/practice_survival_data/dipper.csv")

# We need to modify the Nimble code for the modle in section 3.7
# the observation and transition matrices become time varying, so they receive an additional dimension 
# Also, since there are multiple cohorts of animals marked for the first time each year, 
# the loop must start form the first capture of each individual

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

# We extract detections and non-detections from the data 
y <- dipper %>%
  select(year_1981:year_1987) %>%
  as.matrix()

# We get the occasion of first capture for all individuals 
first <- apply(y, 1, function(x) min(which(x !=0)))

# Now we specify the constants 
my.constants <- list(N = nrow(y),   # number of animals
                     T = ncol(y),   # number of sampling occasions
                     first = first)

# We put the data in a list and add a 1 to the data to code non-detections as 1 and detections as 2
my.data <- list(y = y + 1)

# Now we write a function for the initial values 
# here, for latent states, we assume that individuals are allive through tthe study period 

zinits <- y + 1 # non-detection -> alive
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

# WE can check mixing for the last survivalprobability


