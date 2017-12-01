rm(list=ls())
library(sp)
library(rgdal)
library(spgwr)
sirsam <- readOGR(dsn = "data/geochem_sites.shp")

model.bw <- gwr.sel(Fe_ppm_imp ~ Dlats + DLong , data = sirsam, longlat = TRUE)
model.gauss <- gwr(
  Fe_ppm_imp ~ Dlats + DLong ,
  data = sirsam,
  bandwidth = model.bw,
  se.fit = TRUE,
  longlat = TRUE
)

model.pred <- gwr(
  Fe_ppm_imp ~ Dlats + DLong ,
  data = sirsam,
  bandwidth = model.bw,
  se.fit = TRUE,
  fittedGWRobject = model.gauss,
  longlat = TRUE,
  predictions = TRUE
)

plot(sirsam$Fe_ppm_imp, model.pred$SDF$pred)
