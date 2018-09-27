# A mechanistic model to scale up biophysical processes into geographical size gradients in ectotherms
## Rubalcaba J.G, Gouveia, S.F. & Olalla-Tárraga, M.A.                                                      
##  -> Contact: jg.rubalcaba@gmail.com (Juan G. Rubalcaba)                                                  


 The function theatmodel() is a transient-heat biophysical model that derives the body temperature (Tb) and rate of 
 evaporative water loss  from the skin (EWL) of ectotherms by computing the heat exchanged through: solar radiation, convection, 
 conduction to ground, evaporative cooling, thermal radiation.

 theatmodel() derives body temperature iteratively (through time steps of delta seconds).
 Inputs (outlined below):  morphology/thermal physiology and microcliamtic data. 
 Outputs: Tb after a time step of delta seconds, and the instantaneous rate of EWL.

 Model parameterization:
 Morphology/thermal physiology: we simulate lizard-like and anuran-like ectotherms using allometric functions to
 compute surface area of skin as a function of their body mass, and simulated their thermal physiology using an
 assymetric thermal performance function.

 Microenvironmental data: solar radiation, air and ground temperature, wind speed and relative humidity were computed
 using NicheMapR (Kearney and Porter 2017, Ecography) in the Nearctic and Western Palearctic.
