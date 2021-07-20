#### Script to run power analysis simulation scenario



#### Data simulation
for(reps in 1:nsimdata){
  ## Folders
  dir.create(paste(scen_folder, "/rep", reps, sep = ""))
  
  
  ## Population
  wolverine_pop <- createWolverineTimeSeries(n_0=25,years=10,r=R, stepPars=stepPar,anglePars=anglePar,
                          betaPars=betaPar,stepDist,angleDist,bySex=T,spatial_cov=raster_out,nmoves=psulength*2, movefunc="BVN")
  # This creates a wolverine population with initial population size 25, 10 years of growth with annual rate R, movement governed by step/angle/beta parameters
  # different sexes/classes have different movement patterns, with 2 movements per day. Adults have bivariate normal movement, subadults have correlated random walk
  
  ## Simulate detection
  nc_cameras <- getCameraPoints(poly_out,ceiling(length(poly_out)*prop_cells),0)
  obs_detects <- simDetectByDistance(poly_grid=poly_out,cameraPoints = nc_cameras,animal_occudata = wolverine_pop,p=1,buffer=2000, time_to_switch=NULL)
  # Simulate detection with perfect detection within a 2 km buffer of each camera
  
  ## Prep data for occupancy model
  obs_detects$newtime <- obs_detects$time/2
  # Bin detections
  count <- 1
  dets_to_bin <- numeric(nrow(obs_detects)) 
  while(count <= nrow(obs_detects)-1){
    if(obs_detects$camcell[count]==obs_detects$camcell[count+1] & obs_detects$newtime[count+1]-obs_detects$newtime[count]==0.5){
      dets_to_bin[count+1] <- 1
    }
    count <- count+1
  }
  
  obs_detects <- obs_detects[which(dets_to_bin==0),]
  obs_detects$newtime <- ceiling(obs_detects$newtime)
  
  y <- array(0, dim=c(nrow(nc_cameras), nyears, ntimes))
  for(i in 1:nrow(nc_cameras)){
    for(j in 1:nyears){
      siteyrdets <- obs_detects[which(obs_detects$camcell==nc_cameras$cellID[i] & obs_detects$year==j),]
      detsbyusetime <- cut(siteyrdets$newtime, breaks = seq(0, psulength, by=ssulength))
      y[i,j,] <- as.numeric(table(detsbyusetime))
    }
  }
  
  
  # Set up L - nsites by nyears by ntimes, but all same value
  L <- array(ssulength, dim=c(nrow(nc_cameras), nyears, ntimes))
  
  
  ## Fit occupancy models
  # setup for both models
  win.data <-list(y=y, L=L, nsites=nrow(nc_cameras), nyears=nyears,ntimes=ntimes)
  zinits <- ifelse(apply(y,c(1,2),sum)>0, 1, 0)
  uinits <- ifelse(y > 0, 1, 0)
  
  # implicit dynamic model
  inits=function(){list(int_psi=runif(nyears),int_q=runif(nyears),int_lambda=runif(nyears),z=zinits,u=uinits)}
  params <- c("int_psi", "int_q", "int_lambda", "psiq")
  ni <- 5000; nt <- 2; nb <- 2000; nc <- 3 # 5000;2;2000;3
  out_pwr_impl <- jags(win.data, inits = inits, params, f, n.chains=nc, n.iter=ni, n.burn=nb, n.thin=nt, parallel = TRUE, n.cores = 3)
  
  # explicit dynamic model
  inits=function(){list(int_psi=runif(1),eps=runif(1),gamma=runif(1),int_q=runif(nyears),int_lambda=runif(nyears),z=zinits,u=uinits)}
  params <- c("int_psi", "eps", "gamma", "int_q", "int_lambda", "psi_est", "psiq")
  ni <- 6000; nt <- 2; nb <- 2000; nc <- 3 # 5000;2;2000;3
  out_pwr_expl <- jags(win.data, inits = inits, params, g, n.chains=nc, n.iter=ni, n.burn=nb, n.thin=nt, parallel = TRUE, n.cores = 3)
  
  save(out_pwr_impl, file = paste(scen_folder, "/rep", reps, "/pwrresultsimpl.Rdata", sep = ""))
  save(out_pwr_expl, file = paste(scen_folder, "/rep", reps, "/pwrresultsexpl.Rdata", sep = ""))
}