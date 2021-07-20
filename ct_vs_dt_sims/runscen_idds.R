#library(jagsUI)
###### Small simulation

# Simulate and break up data for the 2 different models
ctdo_data <- list() # continuous-time multi-scale
do_data <- list() # discrete-time multi-scale (3 survey occasions)
for(reps in 1:nsimdata){
  temp_data <- dynocc_PP_sim(nsites = nsites, ntimes = ntimes, psulength = psulength, lambda1 = lam, useprob=useprob, occprob = occprob) # Simulate data
  dir.create(paste(scen_folder, "/rep", reps, sep = ""))
  save(temp_data, file = paste(scen_folder, "/rep", reps, "/simdata.Rdata", sep = ""))
  ctdo_data[[length(ctdo_data)+1]] <- temp_data$ctdo_in
  do_data[[length(do_data)+1]] <- temp_data$do_in
}

### Fit CTDO models

ctdo_results <- list()

#setwd("C:/Users/Robbie/Desktop/wolverines/wolverinepractice/old_and_backup")

# Don't need to concatenate files here? Or do you? No yeah you don't, just need the new data!

for(reps in 1:nsimdata){
  win.data.ctdo <- ctdo_data[[reps]]
  
  uinit <- ifelse(win.data.ctdo$y > 0, 1, 0) # initial use values for JAGS
  zinit <- ifelse(rowSums(win.data.ctdo$y)>0,1,0) # initial occupancy values for JAGS
  
  inits <- function(){list(int_psi=runif(1),int_lambda=runif(1,0,1),int_q=runif(1),u=uinit,z=zinit)}
  
  params <- c("int_psi", "int_lambda", "int_q", "psiq")
  
  ni <- 5000; nt <- 2; nb <- 2000; nc <- 3
  print(paste("Simulating ", reps))
  out_multiyear <- jags(win.data.ctdo, inits = inits, params, ctjags, n.chains=nc, n.iter=ni, n.burn=nb, n.thin=nt, parallel = TRUE, n.cores = 3)
  ctdo_results[[length(ctdo_results)+1]] <- out_multiyear
  # save
  save(out_multiyear, file = paste(scen_folder, "/rep", reps, "/ctdoresults.Rdata", sep = ""))
}

#out_temp <- ctdo_results[[1]]


#########################################################
#########################################################
### Fit do models

do_results <- list()

for(reps in 1:nsimdata){
  win.data.do <- do_data[[reps]]
  
  uinit <- ifelse(apply(win.data.do$x,c(1,2),sum) > 0, 1, 0) # Initial use values for JAGS
  # If any survey occasion within a primary or secondary occasion has a detection, set to 1, else 0
  zinit <- ifelse(apply(win.data.do$x, 1, sum)>0,1,0)

  inits <- function(){list(int_psi=runif(1),int_p=runif(1,0,1),int_q=runif(1),u=uinit,z=zinit)}
  
  params <- c("int_psi", "int_p", "int_q", "psiq") # Note: I added this psiq in after running my simulations the first time. Hopefully won't need to re-run them!
  
  ni <- 5000; nt <- 2; nb <- 2000; nc <- 3
  print(paste("Simulating ", reps))
  out_multiyear <- jags(win.data.do, inits = inits, params, dtjags, n.chains=nc, n.iter=ni, n.burn=nb, n.thin=nt, parallel = TRUE, n.cores = 3)
  do_results[[length(do_results)+1]] <- out_multiyear
  save(out_multiyear, file = paste(scen_folder, "/rep", reps, "/doresults.Rdata", sep = ""))
}