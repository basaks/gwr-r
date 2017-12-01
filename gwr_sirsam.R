rm(list=ls())
library(sp)
library(rgdal)
library(spgwr)
library(parallel)

source('R/read_covariates.R')
sirsam <- readOGR(dsn = "data/geochem_sites.shp")

no.cores <- detectCores() - 1
cl <- makeCluster(no.cores)


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
  fit.points = SpatialPoints(cbind(sirsam$DLong, sirsam$Dlats)),
  cl = cl
)

plot(sirsam$Fe_ppm_imp, model.pred$SDF$pred)
