### "A transient-heat biophyisical model to derive body temperature and rates of evaporative water loss of terrestrial ectotherms" ###
## Juan G. Rubalcaba (Contact: jg.rubalcaba@gmail.com)                                                                              ##
######################################################################################################################################

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
