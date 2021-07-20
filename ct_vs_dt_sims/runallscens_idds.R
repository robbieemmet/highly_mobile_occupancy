#### Script to run simulation study for "Developing an occupancy model for monitoring rare and highly mobile species"
library(jagsUI)
## Change the location of where you want this folder on your computer
suppinfo_local <- "C:/Users//Desktop/wolverines/highlymobile_occ"


##### Run once at beginning of simulation - create data and simulation results folders
# Make sure the tempfolder directory matches where you have supplementary code, then relative path /ct_vs_dt_sims
tempfolder <- paste0(suppinfo_local, "/ct_vs_dt_sims")
for(detectprob in c(0.2, 0.5,0.8)){
  dir.create(paste(tempfolder, "/p", detectprob, sep = ""))
  for(occ in c(0.2, 0.5,0.8)){
    dir.create(paste(tempfolder, "/p", detectprob, "/psi", occ, sep = ""))
    for(q in c(0.2, 0.5,0.8)){
      dir.create(paste(tempfolder, "/p", detectprob, "/psi", occ, "/q", q, sep = ""))
    }
  }
}

# Write continuous-time multi-scale JAGS model file
cat(file={ctjags<-tempfile()}, 
    "
    model{
    
    #Priors
    int_psi ~ dunif(0,1)
    int_q ~ dunif(0,1)
    int_lambda ~ dunif(0,1)
    
    
    #Likelihood
    for(i in 1:nsites){
    psi[i] <- int_psi
    z[i] ~ dbin(psi[i],1)
    for(b in 1:nbaittimes){
    q[i,b] <- int_q
    u[i,b] ~ dbin(z[i]*q[i,b],1)
    y[i,b] ~ dpois(u[i,b]*int_lambda*timebaited[i,b])
    }
    }
    psiq <- int_psi*int_q
    }
    ")

# Write discrete-time multi-scale JAGS model file
cat(file={dtjags<-tempfile()}, 
    "
    model{
    
    #Priors
    int_psi ~ dunif(0,1)
    int_q ~ dunif(0,1)
    int_p ~ dunif(0,1)
    
    
    #Likelihood
    for(i in 1:nsites){
    psi[i] <- int_psi
    z[i] ~ dbin(psi[i],1)
    for(b in 1:nbaittimes){
    q[i,b] <- int_q
    u[i,b] ~ dbin(z[i]*q[i,b],1)
    for(k in 1:3){
    x[i,b,k] ~ dbin(u[i,b]*int_p, 1)
    }
    
    }
    }
    psiq <- int_psi*int_q
    }
    ")

## Function to simulate data
dynocc_PP_sim <- function(nsites = 200, ntimes = 4, nyears = 1, psulength = 120, ssulength = psulength/ntimes, lambda1 = lambda, useprob, occprob){
  lambda <- lambda1
  
  # simulate y - counts at sites
  y <- array(0, dim=c(nsites, ntimes)) # counts at sites
  z_i <- numeric(nsites) # latent occupancy at sites
  u_ib <- array(0, dim=c(nsites, ntimes)) # latent use at site i during time b
  for(i in 1:nsites){
    z_i[i] <- rbinom(1, size = 1, p = occprob)
    for(b in 1:ntimes){
      u_ib[i,b] <- rbinom(1, size = 1, p = useprob)
      # if(sum(u_ib[i,])==0){
      #   z_i[i] <- 0
      # }
      y[i,b] <- rpois(1, z_i[i]*u_ib[i,b]*lambda*ssulength)
    }
  }
  
  L <- matrix(ssulength, nrow=nsites, ncol=ntimes)
  
  # List of CTDO data as output
  ctdo_in <- list(nsites = nsites, nbaittimes = ntimes, y=y, timebaited=L)
  # ^ Add 0 at the beginning then test this!!!
  
  # discrete-time detection history
  nsubreps <- 3 # number of "survey occasions", K in the paper
  subreplength <- ssulength/nsubreps
  max_detects <- max(y) # Maximum number of detections in the data
  # max_detects is used to fill parts of the time array t_ibk with NA's
  t_ibk <- array(NA, dim=c(nsites, ntimes, max_detects)) # times of detections, to be sorted into survey occasions
  x <- array(0, dim=c(nsites, ntimes, nsubreps))
  for(i in 1:nsites){
    for(b in 1:ntimes){
      if(y[i,b]>0){
        temp_t <- runif(y[i,b],0,ssulength) # Detection is a Poisson process so detections are distributed uniformly in time
        t_ibk[i,b,] <- c(temp_t, rep(NA, max_detects-length(temp_t)))
        temp <- cut(t_ibk[i,b,which(!is.na(t_ibk[i,b,]))], breaks = c(0,subreplength,2*subreplength,ssulength))
        temp2 <- temp # note: there's a chance the cutting, if a time of exactly 0 was generated, could have not sorted that into a category and resulted in one continuous-time detection and no discrete-time detections. That never happened in our simulations. To fix, change lowest break in line above to a negative number.
        x[i,b,] <- ifelse(table(temp2) > 0, 1, 0)
      }
    }
  }
  
  
  
  # DO data
  do_in <- list(nsites = nsites, nbaittimes = ntimes, x = x)
  
  
  # return data for continuous-time and discrete-time multi-scale models
  return(list(ctdo_in = ctdo_in, do_in = do_in))
}

## Simulation setup and run
nsites <- 200
psulength <- 120 # Length of primary sampling occasion
ntimes <- 4 # Number of secondary occasions
ssulength <- psulength/ntimes # Length of secondary occasions
nyears <- 1
K <- 3 # Number of tertiary occasions for discrete-time model
nsimdata <- 2 # Number of data sets to simulate per scenario - original is 200
ncores <- 3 # Note: change this if you have fewer than 3 cores available
set.seed(34)
st <- proc.time()[3]
for(detectprob in c(0.2,0.5, 0.8)){
  for(occ in c(0.2, 0.5,0.8)){
    for(q in c(0.2, 0.5,0.8)){
      p <- detectprob # varies
      occprob <- occ # varies
      useprob <- q # varies
      lam <- (-K/ssulength)*log(1-p) # Convert detection probability to rate for Poisson models
      # Specific scenario folder to put rep folders in
      scen_folder <- paste(tempfolder, "/p", detectprob, "/psi", occ, "/q", q, sep = "")
      source(paste(tempfolder, "/runscen_idds.R", sep = ""))
    }
  }
}
et <- proc.time()[3]

et - st

