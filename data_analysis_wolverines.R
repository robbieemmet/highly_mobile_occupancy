#### Data analysis script for wolverines ####
## Updated December 10th, 2020

library(stringr)
library(lubridate)
library(jagsUI)
library(dplyr)

data_folder <- "C:/Users/Robbie/Desktop/wolverines/data/RforRobbie/RforRobbie/"
# ^ Change to where Rdata files are

## Load data
load(paste(data_folder,"detections_Long.Rdata",sep=""))
load(paste(data_folder,"detections_Long.Rdata",sep=""))

ntimes <- 4

time_scaler <- 1

# Add ungroup - 11/30/2020
cams.processed <- cams.processed %>% ungroup()
dets.processed <- dets.processed %>% ungroup()

# Reformat dates
cams.processed$ActiveStart <- as.POSIXct(cams.processed$ActiveStart, format = "%m/%d/%Y %H:%M:%S", tz = "US/Pacific")
cams.processed$ActiveEnd <- as.POSIXct(cams.processed$ActiveEnd, format = "%m/%d/%Y %H:%M:%S", tz = "US/Pacific")
dets.processed$ActiveStart <- as.POSIXct(dets.processed$ActiveStart, format = "%m/%d/%Y %H:%M:%S", tz = "US/Pacific")
dets.processed$ActiveEnd <- as.POSIXct(dets.processed$ActiveEnd, format = "%m/%d/%Y %H:%M:%S", tz = "US/Pacific")
dets.processed$ImageDate <- as.POSIXct(dets.processed$ImageDate, format = "%m/%d/%Y %H:%M:%S", tz = "US/Pacific")

# Create copies of ActiveStart and ActiveEnd for bait model data processing
dets.processed$ActiveStartOld <- dets.processed$ActiveStart
dets.processed$ActiveEndOld <- dets.processed$ActiveEnd
cams.processed$ActiveStartOld <- cams.processed$ActiveStart
cams.processed$ActiveEndOld <- cams.processed$ActiveEnd

# Get rid of camera rows that have no ActiveStart or ActiveEnd
empty_times <- which(is.na(cams.processed$ActiveStart))
cams.processed <- cams.processed[-empty_times,]

## WA only subset
wa.cams <- subset(cams.processed, state=="WA")
wa.dets <- subset(dets.processed, state=="WA")

## Set up use periods - each month when sampling occurred?

cuttimes <- seq(ISOdate(2016, 12, 1, hour = 0, tz = "US/Pacific"), length.out = ntimes+1, by = "30 days")


cutseas <- 1:4

# Some auxiliary variables defining dimensions
nsites <- length(unique(wa.cams$LocationName))
ssulength <- 30
psulength <- ssulength*ntimes

# Create L - effort for each camera

L <- matrix(0, nrow=nsites, ncol=ntimes)
for(i in 1:nsites){
  # Check time zones!
  sitecams <- subset(wa.cams, LocationName==unique(wa.cams$LocationName)[i])
  if(max(sitecams$start) < min(cuttimes) & min(sitecams$end) > max(cuttimes)){
    L[i,] <- rep(ssulength, ntimes)
  } else {
    # this is if any part of the x-month survey period isn't sampled
    
    # beginning time
    L[i,1] <- ssulength-(max(sitecams$start) - min(cuttimes))
    # middle times
    L[i, 2:(ntimes-1)] <- ssulength
    # end time
    L[i, ntimes] <- min(as.numeric(ssulength-(max(cuttimes)-min(sitecams$end))), ssulength)
  }
}


## JAGS models to run the analyses below

ynew <- matrix(0, nrow=nsites, ncol=ntimes)
for(i in 1:nsites){
  sitecams <- subset(wa.cams, LocationName==unique(wa.cams$LocationName)[i])
  sitedets <- subset(wa.dets, LocationName==unique(wa.cams$LocationName)[i])
  if(nrow(sitedets) > 0){
    newtimes <- as.numeric(difftime(sitedets$ImageDate,min(cuttimes),units="mins"))
    
    sitedets$event <- c(1, rep(NA, nrow(sitedets)-1))
    count <- 2
    event <- 1
    while(count <= nrow(sitedets)){
      # Go through newtimes, track time somehow
      current_time <- newtimes[count] # current time
      ongoing_event <- min(which(sitedets$event==event)) # the first time in the current event
      timediff <- current_time-newtimes[ongoing_event] # calculate difference between this time and beginning of event
      if(timediff>60){
        event <- event+1 # reset event
      }
      sitedets$event[count] <- event
      count <- count+1
      
    }
    
    sitedets_unique <- sitedets %>% distinct(event,.keep_all=TRUE)
    occasioncuts <- cut(sitedets_unique$ImageDate, breaks=cuttimes)
    ynew[i,] <- table(occasioncuts)
    
  }
}


## Continuous-time model
# setwd
setwd("C:/Users/Robbie/Desktop/wolverines/wolverinepractice/old_and_backup")
# ^ Set this to where JAGS files are
## impl
win.data <-list(y=ynew, L=L, nsites=nsites,ntimes=ntimes)
zinits <- ifelse(apply(ynew,1,sum)>0, 1, 0)
uinits <- ifelse(ynew > 0, 1, 0)
inits=function(){list(int_psi=runif(1),int_q=runif(1),int_lambda=runif(1),z=zinits,u=uinits)}
params <- c("int_psi", "int_q", "int_p", "psiq")
ni <- 50000; nt <- 10; nb <- 20000; nc <- 3 # 5000;2;2000;3
out_wa_ctdo <- jags(win.data, inits = inits, params, "wolverinepoiswa.txt", n.chains=nc, n.iter=ni, n.burn=nb, n.thin=nt, parallel = TRUE, n.cores = 3)
# Should up the thin rate maybe - done
#print(out_wa_ctdo)

#### Discrete-time data processing (use the below)
K <- 3
xnew <- array(0, dim=c(nsites, ntimes, K))
for(i in 1:nsites){
  sitecams <- subset(wa.cams, LocationName==unique(wa.cams$LocationName)[i])
  sitedets <- subset(wa.dets, LocationName==unique(wa.cams$LocationName)[i])
  if(nrow(sitedets) > 0){
    
    # Start with the beginning of the primary occasion
    newtimes <- as.numeric(difftime(sitedets$ImageDate,min(cuttimes),units="mins"))
    sitedets$event <- c(1, rep(NA, nrow(sitedets)-1))
    count <- 2
    event <- 1
    while(count <= nrow(sitedets)){
      # Go through newtimes, track time somehow
      current_time <- newtimes[count] # current time
      ongoing_event <- min(which(sitedets$event==event)) # the first time in the current event
      timediff <- current_time-newtimes[ongoing_event] # calculate difference between this time and beginning of event
      if(timediff>60){
        event <- event+1 # reset event
      }
      sitedets$event[count] <- event
      count <- count+1
      
    }
    
    sitedets_unique <- sitedets %>% distinct(event,.keep_all=TRUE)
    cuttimes2 <- seq(ISOdate(2016, 12, 1, hour = 0, tz = "US/Pacific"), length.out = ntimes*K+1, by = "10 days")
    occasioncuts <- cut(sitedets_unique$ImageDate, breaks=cuttimes2)
    dettab <- table(occasioncuts)
    for(j in 1:K){
      xnew[i,j,] <- ifelse(dettab[(1:K)+(j-1)*K]>0,1,0)
    }
  }
}

win.data <-list(x=xnew, K=K, nsites=nsites,ntimes=ntimes)
zinits <- ifelse(apply(ynew,1,sum)>0, 1, 0)
uinits <- ifelse(ynew > 0, 1, 0)
inits=function(){list(int_psi=runif(1),int_q=runif(1),int_p=runif(1),z=zinits,u=uinits)}
params <- c("int_psi", "int_q", "int_p", "psiq")
ni <- 50000; nt <- 10; nb <- 20000; nc <- 3 # 5000;2;2000;3
out_wa_dtdo <- jags(win.data, inits = inits, params, "wolverinedtdowa.txt", n.chains=nc, n.iter=ni, n.burn=nb, n.thin=nt, parallel = TRUE, n.cores = 3)
# Should up the thin rate maybe - done
#print(out_wa_dtdo)