# Modelling body temperature of ectotherms           
## Juan G. Rubalcaba                                              

This is a transient-heat model that simulates changes in body temperature of terrestrial ectotherms as a function of solar radiation, convection, conduction to ground, evaporative cooling, and thermal radiation.

````RCode.R```` contains the function ````theatmodel()```` which computes both body temperature and evaporative water loss used in Rubalcaba et al. (2019, Glob Ecol Biogeogr), Rubalcaba et al. (2019, Am Nat) and Rubalcaba & Olalla-Tárraga (2020, J Anim Ecol). 

The function ````theatmodelCpp()```` is just a RCpp version of the above, which makes the calculations more quickly.

Any questions of comments, please contact me jg.rubalcaba@gmail.com
