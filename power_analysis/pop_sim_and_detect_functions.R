# Grid-making function make_grid from http://strimas.com/spatial/hexagonal-grids/
# Correlated random walk code adapted from http://wiki.cbr.washington.edu/qerm/index.php/R/Correlated_Random_Walk
# Note: functions are listed in the order in which they are called in a given simulation run

# Install any of the below packages you do not have installed


#library(dplyr)
#library(plyr)
#library(tidyr)
library(sp)
library(raster)
library(rgeos)
#library(rgbif)
#library(viridis)
#library(gridExtra)
#library(rasterVis)
library(ggplot2)
library(momentuHMM)
#library(spsurvey)
#library(rgdal)
library(jagsUI)
#library(circular)
library(CircStats)
#library(spdep)

make_grid <- function(x, cell_diameter, cell_area, clip = FALSE) {
  # this function creates a hexagonal grid of cells over a raster
  if (missing(cell_diameter)) {
    if (missing(cell_area)) {
      stop("Must provide cell_diameter or cell_area")
    } else {
      cell_diameter <- sqrt(2 * cell_area / sqrt(3))
    }
  }
  ext <- as(extent(x) + cell_diameter, "SpatialPolygons")
  proj4string(ext) <- proj4string(x)
  # generate array of hexagon centers
  g <- spsample(ext, type = "hexagonal", cellsize = cell_diameter, 
                offset = c(0.5, 0.5))
  # convert center points to hexagons
  g <- HexPoints2SpatialPolygons(g, dx = cell_diameter)
  # clip to boundary of study area
  if (clip) {
    g <- gIntersection(g, x, byid = TRUE)
  } else {
    g <- g[x, ]
  }
  # clean up feature IDs
  row.names(g) <- as.character(1:length(g))
  return(g)
}

makeCameraGrid <- function(raster_grid, cell_area){
  # Creates cells of a given size (in meters^2)
  e <- extent(raster_grid)
  temp_poly <- Polygon(rbind(c(e[1],e[3]), c(e[1],e[4]), c(e[2],e[4]), c(e[2],e[3])))
  temp_polys <- Polygons(list(temp_poly), "borders")
  temp_space <- SpatialPolygons(list(temp_polys))
  proj4string(temp_space) <- proj4string(raster_grid)
  poly_grid <- make_grid(temp_space, cell_area = cell_area)
  return(poly_grid)
}

simulateWolverines <- function(n_animals, home_ranges, stepPars, anglePars, betaPars, stepDist, angleDist, bySex = T, givenSexes = T, nf = NULL, nm = NULL, ntm = NULL, spatial_cov, nmoves, movefunc = c("HMM", "BVN", "EXP")){
  if(bySex){ # currently only coded for setting TRUE
    if(!givenSexes){ # Note: any numbers of males, females, transient males you give will not be used unless you manuualy set givenSexes = FALSE
      sexes <- sample(c("female", "male", "transient male"), size=n_animals, prob = c(0.4, 0.2, 0.4),replace=T)
      n_female <- length(which(sexes=="female"))
      n_male <- length(which(sexes=="male"))
      n_tmale <- length(which(sexes=="transient male"))
    } else { # if you provide numbers of adult females, adult males, and subadult transient males
      n_female <- nf
      n_male <- nm
      n_tmale <- ntm
    }
    if(n_female > 0){
      fhrs <- home_ranges[which(home_ranges$sex=="female"),]
      females_raw <- data.frame()
      for(i in 1:nrow(fhrs)){
        # for each female home range, simulate an animal there
        if(movefunc=="HMM"){  
          f_to_add <- simData(nbAnimals=1,nbStates=1,dist=list(step=stepDist,angle=angleDist),
                              Par=list(step=stepPars[1,],angle=anglePars[1,]), beta = matrix(betaPars[1,],nrow=1,ncol=0,byrow=T),lambda=NULL,
                              obsPerAnimal=nmoves, states=TRUE, initialPosition = matrix(c(fhrs[i, "x"],fhrs[i,"y"]), ncol=2), centers = matrix(c(fhrs[i, "x"],fhrs[i,"y"]), ncol=2))
        }
        if(movefunc=="BVN"){
          # I assume uncorrelated x and y
          f_to_add <- data.frame(x=rnorm(nmoves, fhrs[i, "x"], stepPars[1,2]), y=rnorm(nmoves, fhrs[i,"y"],stepPars[1,2]))
        }
        f_to_add$ID <- rep(i, nrow(f_to_add))
        females_raw <- rbind(females_raw,f_to_add)
      }
    } else { # if no females are to be simulated
      females_raw <- data.frame()
    }
    
    if(n_male > 0){
      mhrs <- home_ranges[which(home_ranges$sex=="male"),]
      males_raw <- data.frame()
      for(i in 1:nrow(mhrs)){
        if(movefunc=="HMM"){
          m_to_add <- simData(nbAnimals=1,nbStates=1,dist=list(step=stepDist,angle=angleDist),
                              Par=list(step=stepPars[2,],angle=anglePars[2,]), beta = matrix(betaPars[2,],nrow=1,ncol=0,byrow=T), lambda=NULL,
                              obsPerAnimal=nmoves, states=TRUE, initialPosition = matrix(c(mhrs[i, "x"],mhrs[i,"y"]), ncol=2), centers = matrix(c(mhrs[i, "x"],mhrs[i,"y"]), ncol=2))
        }
        if(movefunc=="BVN"){
          # I assume uncorrelated x and y
          m_to_add <- data.frame(x=rnorm(nmoves, mhrs[i, "x"], stepPars[2,2]), y=rnorm(nmoves, mhrs[i,"y"],stepPars[2,2]))
        }
        m_to_add$ID <- rep(i, nrow(m_to_add))
        males_raw <- rbind(males_raw,m_to_add)
      }
    } else { # if no males are to be simulated
      males_raw <- data.frame()
    }
    
    if(n_tmale > 0){
      tmhrs <- home_ranges[which(home_ranges$sex=="transient male"),]
      tmales_raw <- data.frame()
      for(i in 1:nrow(tmhrs)){
        if(movefunc=="HMM"){
          tm_to_add <- simData(nbAnimals=1,nbStates=1,dist=list(step=stepDist,angle=angleDist),
                               Par=list(step=stepPars[3,],angle=anglePars[3,]), beta = matrix(betaPars[3,],nrow=1,ncol=0,byrow=T), nbCovs=0,lambda=NULL,
                               obsPerAnimal=nmoves, states=TRUE, initialPosition = matrix(c(tmhrs[i, "x"],tmhrs[i,"y"]), ncol=2), centers = matrix(c(tmhrs[i, "x"],tmhrs[i,"y"]), ncol=2))
        }
        if(movefunc=="BVN"){
          # I assume uncorrelated x and y
          tm_steps <- rgamma(nmoves, shape = (stepPars[3,1]^2)/(stepPars[3,2]^2), rate = (stepPars[3,1])/(stepPars[3,2]^2))
          tm_angles <- as.numeric(rwrpcauchy(nmoves, location=anglePars[3,1], rho = anglePars[3,2]))
          tot_angles <- cumsum(tm_angles)
          dX <- tm_steps*cos(tot_angles)
          dY <- tm_steps*sin(tot_angles)
          tm_x <- cumsum(dX)+tmhrs[i, "x"]
          tm_y <- cumsum(dY)+tmhrs[i,"y"]
          tm_to_add <- data.frame(x = tm_x, y = tm_y)
        }
        tm_to_add$ID <- rep(i, nrow(tm_to_add))
        tmales_raw <- rbind(tmales_raw,tm_to_add)
      }
    } else { # if no transient males are to be simulated
      tmales_raw <- data.frame()
    }
    
    females_raw$time <- rep(1:nmoves,n_female)
    males_raw$time <- rep(1:nmoves, n_male)
    tmales_raw$time <- rep(1:nmoves, n_tmale)
    
    if(!is.null(females_raw$x)){
      females_raw <- females_raw[, c("ID", "time", "x", "y")]
    } else {
      females_raw <- data.frame()
    }
    if(!is.null(males_raw$x)){
      males_raw <- males_raw[, c("ID", "time", "x", "y")]
    } else {
      males_raw <- data.frame()
    }
    if(!is.null(tmales_raw$x)){
      tmales_raw <- tmales_raw[, c("ID", "time", "x", "y")]
    } else {
      tmales_raw <- data.frame()
    }
    
    all_raw <- rbind(females_raw, males_raw, tmales_raw)
    prepped_data <- all_raw
    prepped_data$sex = rep(c("female", "male", "transient male"), c(nrow(females_raw), nrow(males_raw), nrow(tmales_raw)))
    data_to_place <- prepped_data
    data_to_place$uniqueID <- paste(data_to_place$sex, as.character(data_to_place$ID), sep="")
    data_to_return <- data_to_place
    return(data_to_return)
  }
}

createWolverineTimeSeries <- function(n_0, years, r, stepPars, anglePars, betaPars, stepDist, angleDist, bySex = T, spatial_cov, nmoves, movefunc=c("HMM", "BVN", "EXP")){
  # r is growth rate in standard discrete-time exponential growth model
  landscape_extent <- extent(spatial_cov)
  init_sexes <- sample(c("female", "male", "transient male"), size=n_0, prob = c(0.4, 0.2, 0.4),replace=T)
  home_ranges <- data.frame(cbind(runif(n_0,landscape_extent[1], landscape_extent[2]), runif(n_0,landscape_extent[3], landscape_extent[4]), init_sexes))
  names(home_ranges) <- c("x", "y", "sex")
  home_ranges$x <- as.numeric(as.character(home_ranges$x))
  home_ranges$y <- as.numeric(as.character(home_ranges$y))
  home_ranges$uniqueID <- numeric(nrow(home_ranges))
  
  #extent(spatial_cov) <- extent(spatial_cov) + c(-15000,15000,-15000,15000)
  # ^ This is only needed if you are using simData (movement using hidden Markov models via momentuHMM)
  
  timeseries <- data.frame()
  for(year in 1:years){
    n_current <- floor(n_0*r^year)
    if (n_current == 0){ # population is 0
      return(timeseries)
      break
    }
    if(r <= 1){
      current_hr_length <- nrow(home_ranges) # this the number of home ranges/wolverines from the previous year
      home_ranges <- home_ranges[sort(sample(1:current_hr_length,n_current)),] # pick some of the previous home ranges (n_current of them)
      sexes <- home_ranges[,3]
      fs <- which(sexes=="female")
      ms <- which(sexes=="male")
      tms <- which(sexes=="transient male")
      uniqueids <- numeric(nrow(home_ranges))
      uniqueids[fs] <- 1:length(fs) # Note: these unique IDs change across years
      uniqueids[ms] <- 1:length(ms)
      uniqueids[tms] <- 1:length(tms)
      home_ranges$uniqueID <- paste(home_ranges$sex, uniqueids, sep="")
      # sample new home ranges for transients
      home_ranges[tms, c("x", "y")] <- cbind(runif(length(tms),landscape_extent[1], landscape_extent[2]),
                                             runif(length(tms),landscape_extent[3], landscape_extent[4]))
    } else {
      # What to do if the population is growing (or staying stable)
      current_hr_length <- nrow(home_ranges)
      animals_to_add <- n_current - current_hr_length
      new_sexes <- sample(c("female", "male", "transient male"), size=animals_to_add, replace=T) # Note: new sexes are added uniformly, not in proportion to starting population sex ratios
      new_home_ranges <- data.frame(cbind(runif(animals_to_add,landscape_extent[1], landscape_extent[2]), runif(animals_to_add,landscape_extent[3], landscape_extent[4]), new_sexes))
      names(new_home_ranges) <- c("x", "y", "sex")
      new_home_ranges$x <- as.numeric(as.character(new_home_ranges$x))
      new_home_ranges$y <- as.numeric(as.character(new_home_ranges$y))
      new_home_ranges$uniqueID <- numeric(nrow(new_home_ranges))
      
      home_ranges <- rbind(home_ranges, new_home_ranges)
      sexes <- home_ranges[,3]
      fs <- which(sexes=="female")
      ms <- which(sexes=="male")
      tms <- which(sexes=="transient male")
      uniqueids <- numeric(nrow(home_ranges))
      uniqueids[fs] <- 1:length(fs)
      uniqueids[ms] <- 1:length(ms)
      uniqueids[tms] <- 1:length(tms)
      home_ranges$uniqueID <- paste(home_ranges$sex, uniqueids, sep="")
      home_ranges[tms, c("x", "y")] <- cbind(runif(length(tms),landscape_extent[1], landscape_extent[2]),
                                             runif(length(tms),landscape_extent[3], landscape_extent[4]))
      
      
      
    }
    data_to_add <- simulateWolverines(n_current, home_ranges, stepPars, anglePars, betaPars, stepDist, angleDist, bySex = T, givenSexes = T, nf = length(fs), nm = length(ms), ntm = length(tms), spatial_cov, nmoves = nmoves, movefunc=movefunc)
    data_to_add$year <- rep(year, nrow(data_to_add))
    timeseries <- rbind(timeseries, data_to_add)
  }
  return(timeseries)
}

getCameraPoints <- function(poly_grid, n_sentinels, n_rotating, time_to_switch=NULL){
  # Note: this function uses the centers of cells as cameras, but this could be changed easily
  
  temp <- getSpPPolygonsLabptSlots(poly_grid) # these are the centers of each hexagon, more or less
  camera_sentinels <- sample(c(1:nrow(temp)), n_sentinels) # decide which cameras will be stationary
  camera_rotating <- sample(c(1:nrow(temp))[-camera_sentinels], n_rotating) # decide which cameras will rotate
  camera_cells <- data.frame(rbind(temp[camera_sentinels,], temp[camera_rotating,]))
  camera_cells$cellID <- c(camera_sentinels, camera_rotating)
  names(camera_cells) <- c("x", "y", "cellID")
  return(data.frame(x = camera_cells$x, y=camera_cells$y, camera_class = c(rep("sentinel", n_sentinels), rep("rotating", n_rotating)), cellID = camera_cells$cellID))
}

changeCameraPoints <- function(poly_grid, cameraPoints){
  # Note: this function is only useful if you rotate any cameras
  temp <- getSpPPolygonsLabptSlots(poly_grid)
  # dont_touch tells you which cameras are in the previous year, not to be rotated to
  dont_touch <- which(1:nrow(temp) %in% cameraPoints$cellID)
  if(nrow(temp)-length(dont_touch) > length(which(cameraPoints$camera_class=="rotating"))){
  new_ones <- sample(c(1:nrow(temp))[-dont_touch], length(which(cameraPoints$camera_class=="rotating")))
  cameraPoints$cellID[which(cameraPoints$camera_class=="rotating")] <- new_ones
  cameraPoints[which(cameraPoints$camera_class=="rotating"), c("x", "y")] <- temp[new_ones,]
  } else {
    new_ones <- sample(c(1:nrow(temp))[-dont_touch], nrow(temp)-length(dont_touch))
    rotating_cams <- which(cameraPoints$camera_class=="rotating")
    to_rotate <- rotating_cams[1:length(new_ones)]
    cameraPoints$cellID[to_rotate] <- new_ones
    cameraPoints[to_rotate, c("x", "y")] <- temp[new_ones,]
  }
  
  
  return(cameraPoints)
}

simDetectByDistance <- function(poly_grid=NULL, cameraPoints, animal_occudata, p, buffer,time_to_switch=NULL){
  # Simulates detection by a camera if a wolverine goes within "buffer" distance of it
  detects <- numeric(nrow(animal_occudata))
  camcell <- numeric(nrow(animal_occudata))
  #cells_checked <- c() # This tells you which cells you've checked, if you rotate any cameras
  dist_to_camera <- pointDistance(matrix(c(animal_occudata$x, animal_occudata$y),ncol=2), matrix(c(cameraPoints$x, cameraPoints$y),ncol=2), lonlat = FALSE)
  # ^ This is the distance between each wolverine location and each camera, in a UTM coordinate system (meters)
  for(i in 1:nrow(dist_to_camera)){
    # ^ for each row of the wolverine location data
    which_detected <- which(dist_to_camera[i,] <= buffer) # For a given animal observation, which cameras might detect it
    if(length(which_detected) > 0){
      # if that wolverine location is detected/detectable, pick exactly one camera to detect that wolverine at that time
      if(length(which_detected) > 1){
        choose_one <- sample(which_detected, size = 1, replace = T)
      } else {
        choose_one <- which_detected
      }
      camcell[i] <- cameraPoints[choose_one,"cellID"] # get the cellID column from the cameraPoints data frame (output of getCameraPoints) for the camera that detects the wolverine location
      detects[i] <- rbinom(1,1,p) # Wolverine is detected within buffer with probability p
    } else {
      # If no detections, no camera information and not detected
      camcell[i] <- NA
      detects[i] <- 0
    }
  }
  
  animal_occudata <- data.frame(animal_occudata, detects = detects, camcell = camcell)
  data_to_return <- animal_occudata[which(animal_occudata$detects==1),c("time", "camcell", "year")]
  
  
  # Add rows for cameras where no detections occurred
  #cells_without_detections <- cells_checked[!(cells_checked %in% data_to_return$camcell)]
  #empty_adds <- data.frame(time = rep(NA, length(cells_without_detections)), camcell = cells_without_detections, year = rep(NA, length(cells_without_detections)))
  #data_to_return <- rbind(data_to_return, empty_adds)
  return(data_to_return)
  
}
