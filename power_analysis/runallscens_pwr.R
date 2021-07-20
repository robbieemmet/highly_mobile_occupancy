#### Script to run power analysis for "Developing an occupancy model for monitoring rare and highly mobile species"
## Change the location of where you want this folder on your computer
suppinfo_local <- "C:/Users//Desktop/wolverines/highlymobile_occ"

source(paste0(suppinfo_local, "/power_analysis/pop_sim_and_detect_functions.R"))

raster_out <- raster()
extent(raster_out) <- c(622849, 695033.4, 5343033, 5441739)
# This is the extent of an area of the North Cascades
proj4string(raster_out) <- "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
poly_out <- makeCameraGrid(raster_out, cell_area = 225000000) # cell area is 225 km^2 = 225,000,000 m^2

##### Run once at beginning of simulation
tempfolder <- paste0(suppinfo_local, "/power_analysis")
for(R in c(0.933, 0.972, 1.022, 1.041)){
  dir.create(paste(tempfolder, "/R", R, sep = ""))
  for(prop_cells in c(0.2, 0.5, 0.8)){
    dir.create(paste(tempfolder, "/R", R, "/pcells", prop_cells, sep = ""))
  }
}

#### Write model file
## impl
cat(file={f<-tempfile()},
    "
  model{
      
      #Priors
      # between years
      for(j in 1:nyears){
      int_psi[j] ~ dunif(0,1)
      
      # within years
      
      int_q[j] ~ dunif(0,1)

      int_lambda[j] ~ dunif(0,1)
      }
      
      
      #Likelihood
      for(i in 1:nsites){
        # Occupancy 
        for(j in 1:nyears){
        z[i,j] ~ dbinom(psi[i,j],1)
        psi[i,j] <- int_psi[j]
        for(b in 1:ntimes){
          # Use
          u[i,j,b] ~ dbinom(z[i,j]*q[i,j,b],1)
          q[i,j,b] <- int_q[j]
          # Detection
          y[i,j,b] ~ dpois(u[i,j,b]*int_lambda[j]*L[i,j,b])
        }
        }
      }

      for(j in 1:nyears){
        psiq[j] <- int_psi[j]*int_q[j]
      }
    }    
    "
)

## expl
cat(file={g<-tempfile()},
    "
    
    model{
    
    #Priors
    int_psi ~ dunif(0,1)
    eps ~ dunif(0,1)
    gamma ~ dunif(0,1)
    # between years
    for(j in 1:nyears){
    
    # within years
    
    int_q[j] ~ dunif(0,1)
    
    int_lambda[j] ~ dunif(0,1)
    }
    
    #Likelihood
    for(i in 1:nsites){
    # Occupancy - year 1
    z[i,1] ~ dbinom(pi[i,1],1)
    pi[i,1] <- int_psi
    
    for(b in 1:ntimes){
    # Use
    u[i,1,b] ~ dbinom(z[i,1]*q[i,1,b],1)
    q[i,1,b] <- int_q[1]
    # Detection
    y[i,1,b] ~ dpois(u[i,1,b]*int_lambda[1]*L[i,1,b])
    }
    
    for(j in 2:nyears){
    z[i,j] ~ dbinom(pi[i,j],1)
    pi[i,j] <- z[i,j-1]*(1-eps) + (1-z[i,j-1])*gamma
    for(b in 1:ntimes){
    # Use
    u[i,j,b] ~ dbinom(z[i,j]*q[i,j,b],1)
    q[i,j,b] <- int_q[j]
    # Detection
    y[i,j,b] ~ dpois(u[i,j,b]*int_lambda[j]*L[i,j,b])
    }
    }
    }
    
    # Derived quantity - occupancy in each year
    for(j in 1:nyears){
    psi_est[j] <- sum(pi[,j])/nsites
    psiq[j] <- psi_est[j]*int_q[j]
    }
    }    
    
    "
)

set.seed(34)
psulength <- 120 # length of primary sampling occasion
ntimes <- 4 # number of secondary occasions
nyears <- 10
ssulength <- psulength/ntimes # length of secondary occasions
stepPar <- matrix(c(NA,5000,NA,8000,7200,4200), nrow=3, byrow=T) # mean_1, mean_2, sd_1, sd_2, zeromass_1, zeromass_2
anglePar <- matrix(c(NA,NA,NA,NA,3,0.09),nrow=3, byrow=T) # mean_1, mean_2, concentration_1, concentration_2
# ^ NOTE: If using HMM movement models, all of these must be filled
# If using BVN movement for adults and correlated random walk for transients, only entries 2, 4 (BVN standard deviation), 5, 6 (gamma step parameters) must be filled for stepPar
# If using BVN movement for adults and CRW for transients, only entries 5 and 6 (turning angle parameters) must be filled
stepDist <- "gamma" # for HMM movement model only
angleDist <- "wrpcauchy" # for HMM movement model only
betaPar <- matrix(c(2,2,2), nrow=3, ncol = 1, byrow=T)
nsimdata <- 2 # number of simulations per scenario - original is 100
ncores <- 3 # NOTE: You should change this is you have fewer than 3 cores available on your computer
# But you don't need more than 3 as no model has more than 3 chains
st <- proc.time()[3]
for(R in c(0.933, 0.972, 1.022, 1.041)){
  for(prop_cells in c(0.2, 0.5, 0.8)){
      R <- R # varies
      prop_cells <- prop_cells # varies
      # Specific scenario folder to put rep folders in
      scen_folder <- paste(tempfolder, "/R", R, "/pcells", prop_cells, sep = "")
      source(paste(tempfolder, "/runscen_pwr.R", sep = ""))
  }
}
et <- proc.time()[3]

et - st

