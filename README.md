# A transient-heat biophyisical model to derive body temperature and rates of evaporative water loss of terrestrial ectotherms.           
## Juan G. Rubalcaba (jg.rubalcaba@gmail.com)                                              

RCode.R contains the transient-heat biophyisical model used in the manuscripts: Rubalcaba et al. "A mechanistic model to scale up biophysical processes into geographical size gradients in ectotherms", and Rubalcaba et al. "Upscaling microclimatic conditions into body temperature distribution of ectotherms"

The function theatmodel() is a transient-heat biophysical model that derives body temperature (Tb) and rate of 
evaporative water loss from the skin (EWL) of ectotherms by computing the heat exchanged through solar radiation, convection, 
conduction to ground, evaporative cooling, thermal radiation.

The function theatmodelCpp() is the RCpp version of the transient-heat model that allows more rapid iterative computation.
