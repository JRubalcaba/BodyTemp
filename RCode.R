### "A mechanistic model to scale up biophysical processes into geographical size gradients in ectotherms" ###
## Rubalcaba J.G, Gouveia, S.F. & Olalla-Tárraga, M.A. 
##  -> Contact: jg.rubalcaba@gmail.com (Juan G. Rubalcaba)
##################################################################################################################

## The function theatmodel() is a transient-heat biophysical model that derives the body temperature (Tb) and rate of 
# evaporative water loss  from the skin (EWL) of ectotherms by computing the heat exchanged through: solar radiation, convection, 
# conduction to ground, evaporative cooling, thermal radiation.

## theatmodel() derives body temperature iteratively (through time steps of delta seconds).
# Inputs (outlined below):  morphology/thermal physiology and microcliamtic data. 
# Outputs: Tb after a time step of delta seconds, and the instantaneous rate of EWL.

## Model parameterization:
# Morphology/thermal physiology: we simulate lizard-like and anuran-like ectotherms using allometric functions to
# compute surface area of skin as a function of their body mass, and simulated their thermal physiology using an
# assymetric thermal performance function.

# Microenvironmental data: solar radiation, air and ground temperature, wind speed and relative humidity were computed
# using NicheMapR (Kearney and Porter 2017, Ecography) in the Nearctic and Western Palearctic.


theatmodel <- function(Tb,       # Body temperature at time t (the function gives body temperature at time t+delta) (ºC)
                            A,        # Total surface area of skin (m^2)
                            M,        # Body mass (g)
                            Ta,       # Air temperature (ºC)
                            Tg,       # Ground temperature (ºC)
                            S,        # Solar radiation (W m^-2)
                            v,        # Wind speed (m s^-1)
                            HR,       # Relative humidity (%)
                            skin_humidity = 1, # Skin humidity (0-1)
                            r,        # Total skin resitance to water loss (s m^-1) 
                            posture = 1, # postural adjustment: proportion of dorsal skin surface exposed to the sunlight
                            vent = 1/3,  # proportion of ventral A
                            delta = 60, # Number of time steps (seconds)
                            C = 3.6    # Specific heat capacity (J g-1 ºC-1)
                            )
  {
    for(i in 1:delta){
      
      ## Heat transfer function
      Q <- function(Tb){
        qabs <- A1 * a * S + 0.2 * A1 * S + 0.2 * A1 * S # Absorption of solar radiation
        qevap <-  L * EWL_estimate                     # Evaporative cooling
        qrad <- A1 * sigma * epsilon * (Tb^4 - Ta^4)   # Emission of long-wave radiation
        qconv <- A1 * hc * (Tb - Ta)                   # COnvection
        qcond <- Ag * hg * (Tb - Tg)                   # Conduction to the ground
        qnet <- qabs - qevap - qrad - qconv - qcond    # Net heat flux
        return(qnet)
      }
      
      ## Model parameters
      # Animal geometry and color
      Ag = vent * A   # ventral area (m^2)
      A1 = (1-vent) * A * posture # exposed area (m^2)
      l = sqrt(A)     # caracterstic length (m)
      C = C           # average heat capacity body (J g^-1 ºC^-1)
      epsilon = 0.9   # Emissivity to long-wave radiation
      a = 0.8         # Absorptance to short-wave solar radiation
      
      # Constants
      L = 2257         # latent heat of vaporization of water (J g^-1)
      sigma = 5.670367e-8 # Stefan-Boltzmann constant (W m^-2 ºC^-4)
      
      # Estimation convection coefficient
      rho = 101325 / (287.04 * (Ta + 273))
      nu = -1.1555e-14*(Ta+273)^3 + 9.5728e-11*(Ta+273)^2 + 3.7604e-08*(Ta+273) - 3.4484e-06    # kinematic viscosity (m^2 s^-1)
      kf = 1.5207e-11*(Ta+273)^3 - 4.8574e-08*(Ta+273)^2 + 1.0184e-04*(Ta+273) - 3.9333e-04 # thermal conductivity air (W m^-1 ºC^-1)
      
      Re = l * v / nu   # Reynolds number
      c = 0.2           # Constants relating Re and Nu numbers 
      n = 0.7 
      Nu = 2 + c * Re^n # Nusselt number
      hc = Nu * kf / l  # Convection coefficient (W m-2 ºC-1) 
      
      # Estimation of conduction coefficient
      ksub = 0.027        # thermal conductivity of substrate (W ºC-1)
      ts = 0.025 * (0.001 * M / (pi * 1000))^0.2   # thinkness of surface layer (Kleiber 1972)
      hg = ksub / ts    # Conduction coeficient (W m-2 ºC-1) (Stevenson 1985)
      
      # Estimate evaporative water loss (g s-1)
      HR_prop = HR * 0.01 
      ps_a = exp(77.3450 + 0.0057 * (Ta+273) - 7235 / (Ta+273)) / (Ta+273)^8.2  # Air vapor pressure (Pa)
      rho_a = HR_prop * 2.2 * ps_a / (Ta+273)         # Trasform into density (g m-3)
      
      ps_b = exp(77.3450 + 0.0057 * (Tb+273) - 7235 / (Tb+273)) / (Tb+273)^8.2 # Body vapor density (Pa)
      rho_b = skin_humidity * 2.2 * ps_b / (Tb+273) # Trasform into density (g m-3) 
      
      EWL_estimate = A1 * (rho_b - rho_a) / r  # Evaporative water loss (g s-1)
      
      ## Model iteration
      # Body Tempearature at t + 1" 
      Tb = Tb +  1 / (C * M) * Q(Tb)
    }
  
    # Output: return body temperature at time t + delta (ºC) and instantaneous rate of evaporative water loss (g s-1) 
    return(c(Tb, EWL_estimate)) 
    
  }

# Example:

# This is an example of how to use the function "theatmodel()". You can just run the function's script (above) and then run
# the "example" code below to get the outcome: body temperature and rate of evaporative water loss after a time step of delta seconds.

Tb = 20 # starting Tb
M = 10 # Body Mass (g)
A = 11 * M^(2/3) * 1e-4  # Total surface area of skin (m-2)
r = 60000 # Skin resistance to water loss (s m^-1)
C = 3.6  # Specific heat capacity of the body (J g-1)

Ta = 25 # Air temperature (ºC)
Tg = 30 # Ground temperature (ºC)
S = 400 # Total solar radiation (W m-2)
v = 10  # wind speed (m s-1)
HR = 40 # relative humidity (%)

delta = 60 # Number of steps (seconds)

# The function gives two outputs: Tb (ºC) at t + delta, and instantaneous rate of EWL at t + delta (g s-1)
output <- theatmodel(Tb, A, M, Ta, Tg, S, v, HR, skin_humidity = 1, r, posture = 1, vent = 1/3, delta, C)
output[1] # Tb (ºC) after delta seconds
output[2] * 1000 * 3600 / M # instantaneous EWL (mg g-1 h-1) at time t + delta

############################################################
### PARAMETERIZATION I: Animal morphology and physiology ###
############################################################

#### Morphology

min_M <- 5; max_M <- 200 # Range of body masses (g)
n_spec <- 3 # Number of size clases
M = seq(min_M, max_M, length.out = n_spec); S = 11 * M^(2/3) * 1e-4 # Calculate surfaces (m2) and masses (g) with allomatric function
geom <- data.frame(M,S); rm(S, M)

plot(S ~ M, geom, type = "l", ylab = expression(paste("Surface (", m^{2}, ")")), xlab = "Body mass (g)")

#### Thermal physiology

Topt = 33 # Optimal temperature (ºC)
performance <- function(x,Topt) if((x-Topt)<0){exp(-((x-Topt)/(4*2))^2)} else {1-((x-Topt)/(Topt-1*25))^2}
w <- sapply(seq(0,45,0.1), function(x)performance(x,Topt))
plot(w ~ seq(0,45,0.1), type="l", lwd=1, ylim=c(0,1), ylab="Relative performance", xlab="Body Temperature (ºC)")
abline(v=Topt, lty=3)

lambda = 20  # Thermoregulatory constraint
r = 60000 # Total resitance to water loss (s m^-1) -> Reptiles (60000 sm-1; amphibians 300 sm-1)


#################################################
#### PARAMETERIZATION II: Environmental data  ###
#################################################

# Extract microclimatic variables using NicheMapR (Kearney and Porter 2017, Ecography) 

require("raster")

# First, you need to upload a raster with the resolution and extent that you want for your maps
# here, we use 1ºx1º grid cell systems in North America and Western Palearctic

shapeTotal <- raster("....")	# load a raster to build the maps (this is just the geographical support to build the maps, we are not extracting climate data yet)
extent <- c(-130, -55, 20, 70) # Nearctic: c(-130, -55, 20, 70) # Western Palearctic: c(-20, 50, 35, 71) 
map <- crop(shapeTotal, extent)
map <- aggregate(map, fact=1/res(map)[1]) # Adjust spatial resolution 

regions <- shapefile("...") # load shapefiles of biogeographical regions
nearctic <- regions[regions[[5]]=="Nearctic",]

map <- mask(map, nearctic) # Create the raster

ext <- extract(map, nearctic, cellnumbers=TRUE)[[1]][,1]
xy.values <- xyFromCell(map, ext) # This is the dataset containing xy coords of our cells

# Now, we will use NicheMapR to extract microclimatic variables at each cell
# We run NicheMapR 3 times to derive microclimatic variables (solar radiation, air and ground temperature, wind speed,
# and relative humidity) in three microenvironments (full sun, 50% shade and 90% shade)

require(NicheMapR)

month = 6

### Map full sun
data_map_sun <- list()
for(i in 1:length(xy.values[,1])){
  tryCatch(
    micro <- micro_global(loc = xy.values[i,], minshade = 0, maxshade = 100),
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  metout <- as.data.frame(micro$metout) 
  soil <- as.data.frame(micro$soil) 
  
  day = unique(metout$DOY)[month] 
  TIME <- subset(metout, DOY == day)$TIME
  RAD <- subset(metout, DOY == day)$SOLR 
  AIRT <- subset(metout, DOY == day)$TALOC
  SOILT <- subset(soil, DOY == day)$D0cm
  RH <- subset(metout, DOY == day)$RHLOC
  V <- subset(metout, DOY == day)$VLOC
  
  data_map_sun[[i]] <- data.frame(TIME, RAD, AIRT, SOILT, RH, V)
} 
# Change time resolution (minutes)
data_map_sun_min <- list()
modif <- array(NA, dim=c(1440, 6)); modif[,1] <- 1:1440
for(i in 1:length(xy.values[,1])){
  for(j in 2:6){
    f <- approxfun(data_map_sun[[i]][,1], data_map_sun[[i]][,j])
    S <- f(seq(0, 1380, length  = 1440))
    modif[,j] <- S
  }
  data_map_sun_min[[i]] <- modif
}

### Map 50% shade
data_map_mid <- list()
for(i in 1:length(xy.values[,1])){
  tryCatch(
    micro <- micro_global(loc = xy.values[i,], minshade = 50, maxshade = 100),
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  metout <- as.data.frame(micro$metout) 
  soil <- as.data.frame(micro$soil) 
  
  day = unique(metout$DOY)[month] 
  TIME <- subset(metout, DOY == day)$TIME
  RAD <- subset(metout, DOY == day)$SOLR * 0.5
  AIRT <- subset(metout, DOY == day)$TALOC
  SOILT <- subset(soil, DOY == day)$D0cm
  RH <- subset(metout, DOY == day)$RHLOC
  V <- subset(metout, DOY == day)$VLOC
  
  data_map_mid[[i]] <- data.frame(TIME, RAD, AIRT, SOILT, RH, V)
}
# Change time resolution
data_map_mid_min <- list()
modif <- array(NA, dim=c(1440, 6)); modif[,1] <- 1:1440
for(i in 1:length(xy.values[,1])){
  for(j in 2:6){
    f <- approxfun(data_map_mid[[i]][,1], data_map_mid[[i]][,j])
    S <- f(seq(0, 1380, length  = 1440))
    modif[,j] <- S
  }
  data_map_mid_min[[i]] <- modif
}

### Map 90% shade
data_map_shade <- list()
for(i in 1:length(xy.values[,1])){
  tryCatch(
    micro <- micro_global(loc = xy.values[i,], minshade = 90, maxshade = 100),
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  metout <- as.data.frame(micro$metout) 
  soil <- as.data.frame(micro$soil) 
  
  day = unique(metout$DOY)[month] 
  TIME <- subset(metout, DOY == day)$TIME
  RAD <- subset(metout, DOY == day)$SOLR * 0.9
  AIRT <- subset(metout, DOY == day)$TALOC
  SOILT <- subset(soil, DOY == day)$D0cm
  RH <- subset(metout, DOY == day)$RHLOC
  V <- subset(metout, DOY == day)$VLOC
  
  data_map_shade[[i]] <- data.frame(TIME, RAD, AIRT, SOILT, RH, V)
} 
# Change time resolution
data_map_shade_min <- list()
modif <- array(NA, dim=c(1440, 6)); modif[,1] <- 1:1440
for(i in 1:length(xy.values[,1])){
  for(j in 2:6){
    f <- approxfun(data_map_shade[[i]][,1], data_map_shade[[i]][,j])
    S <- f(seq(0, 1380, length  = 1440))
    modif[,j] <- S
  }
  data_map_shade_min[[i]] <- modif
}

## The above code creates 3 lists (sun, 50% shade and 90% shade), which contain one matrix for each cell of the map.  
# Each matrix has 1440 rows (one per minute) and 6 columns: (1) minute of the day, (2) solar rad, (3) Ta, (4) Tg, (5) RH, (6) wind speed

head(data_map_mid_min[[50]])

plot(data_map_sun_min[[50]][,3], ylab = "Air Temperature, ºC", xlab = "Time, min")

## Check the maps
j=0
for(i in ext){
  j = j+1
  map[i] <- max(data_map_sun_min[[j]][,3]) # e.g., max air temperature
}
plot(mask(map, nearctic)); lines(nearctic)



######################
#### RUN THE MODEL ###
######################

# This script uses the transient-heat biophysical model to derive daytime Tb and EWL of 
# thermoregulating ectotherms that select microenvironments within the repertoire.

iter = 1440 # Number of iterations
Tb_map <- list(); EWL_map <- list()
Tb <- array(Topt, dim=c(iter, n_spec)); EWL <- array(NA, dim=c(iter, n_spec))
cells <- length(xy.values[,1]) 
for(k in 1:cells){
  for(i in 1:n_spec){
    A = geom$S[i]
    M = geom$M[i]
    
    for(j in 2:iter){
      S_sun <- data_map_sun_min[[k]][j,2] 
      S_mid <- data_map_mid_min[[k]][j,2] 
      S_shade <- data_map_shade_min[[k]][j,2] 
      
      Ta_sun <- data_map_sun_min[[k]][j,3]
      Ta_mid <- data_map_mid_min[[k]][j,3]
      Ta_shade <- data_map_shade_min[[k]][j,3]
      
      Tg_sun <- data_map_sun_min[[k]][j,4]
      Tg_mid <- data_map_mid_min[[k]][j,4]
      Tg_shade <- data_map_shade_min[[k]][j,4]
      
      HR_sun <- data_map_sun_min[[k]][j,5] 
      HR_mid <- data_map_mid_min[[k]][j,5] 
      HR_shade <- data_map_shade_min[[k]][j,5] 
      
      V_sun <- data_map_sun_min[[k]][j,6] 
      V_mid <- data_map_mid_min[[k]][j,6] 
      V_shade <- data_map_shade_min[[k]][j,6]
      
      # Body temperature and EWL for each microenvironment and posture
      t1_sun <- theatmodel(Tb[j-1,i], A, M, Ta_sun, Tg_sun, S_sun, V_sun, HR_sun, skin_humidity = 1, r, posture = 1, delta = 60)
      t1_mid <- theatmodel(Tb[j-1,i], A, M, Ta_mid, Tg_mid, S_mid, V_mid, HR_mid, skin_humidity = 1, r, posture = 1, delta = 60) 
      t1_shade <- theatmodel(Tb[j-1,i], A, M, Ta_shade, Tg_shade, S_shade, V_shade, HR_shade, skin_humidity = 1, r, posture = 1, delta = 60)
      
      t1_sun_posture <- theatmodel(Tb[j-1,i], A, M, Ta_sun, Tg_sun, S_sun, V_sun, HR_sun, skin_humidity = 1, r, posture = 0.1, delta = 60)
      t1_mid_posture <- theatmodel(Tb[j-1,i], A, M, Ta_mid, Tg_mid, S_mid, V_mid, HR_mid, skin_humidity = 1, r, posture = 0.1, delta = 60)
      t1_shade_posture <- theatmodel(Tb[j-1,i], A, M, Ta_shade, Tg_shade, S_shade, V_shade, HR_shade, skin_humidity = 1, r, posture = 0.1, delta = 60)
      
      # Behavioural thermoregulation
      T_t1 <- rbind(t1_sun[1], t1_mid[1], t1_shade[1],
                    t1_sun_posture[1], t1_mid_posture[1], t1_shade_posture[1])
      
      EWL_t1 <- rbind(t1_sun[2], t1_mid[2], t1_shade[2],
                      t1_sun_posture[2], t1_mid_posture[2], t1_shade_posture[2])
      
      w_i <- sapply(T_t1, function(x) performance(x, Topt)) # Calculates performance for each temperature
      
      p <- exp(lambda*w_i)/sum(exp(lambda*w_i)) # Calculate multinomial probability and accept one microenv with probability p
      select <- rmultinom(1, 1, p)
      
      # Set body temperature and EWL
      
      Tb[j,i] <- T_t1[which(select==1)]
      EWL[j,i] <- EWL_t1[which(select==1)]
    }
  }
  Tb_map[[k]] <- Tb
  EWL_map[[k]] <- EWL
  print(paste0(round(k / cells * 100, 2), " %"))
}

# Now we can build the maps of thermoregulatory accuracy and evaporative water loss of ectotherms 
## Accuracy map

bsize <- 1 # Select body size: 1 (5g), 2 (100g), 3(200g) 
activity_period <- 420:1080 # 7 to 18 h

j=0
for(i in ext){
  j = j+1
  accuracy <- 1 / mean(abs(Tb_map[[j]][activity_period, bsize] - Topt))
  map[i] <- accuracy
}
plot(map)

## EWL map

bsize <- 1 # Select body size: 1 (5g), 2 (100g), 3(200g) 
activity_period <- 420:1080 # 7 to 18 h
j=0
for(i in ext){
  j = j+1
  ewl <- mean(EWL_map[[j]][activity_period, bsize]) 
  map[i] <- ewl 
}
plot(mask(map, nearctic)); lines(nearctic)

