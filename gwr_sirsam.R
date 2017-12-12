rm(list=ls())
library(sp)
library(rgdal)
library(spgwr)
library(parallel)
library(raster)


source('R/read_covariates.R')
sirsam <- readOGR(dsn = "data/geochem_sites.shp")
covariates.file <- 'data/sirsam_covariates_Na.txt'

cl <- makeCluster(detectCores() - 1)

model.bw <- gwr.sel(Na_ppm_imp ~ Dlats + DLong , data = sirsam, longlat = TRUE)

model.gauss <- gwr(
  Na_ppm_imp ~ Dlats + DLong ,
  data = sirsam,
  bandwidth = model.bw,
  se.fit = TRUE,
  longlat = TRUE
)

model.pred <- gwr(
  Na_ppm_imp ~ Dlats + DLong,
  data = sirsam,
  bandwidth = model.bw,
  se.fit = TRUE,
  fittedGWRobject = model.gauss,
  longlat = TRUE,
  fit.points = SpatialPoints(cbind(sirsam$DLong, sirsam$Dlats)),
  cl = cl
)

# basic regression prediction plot using only lat lon
plot(sirsam$Na_ppm_imp, model.pred$SDF$pred)

# pick a sample raster for gathering the indices of the targets in rasters
covariates.list = ParseText(covariates.file)
sample.ras <- raster(covariates.list[1])
sample.crs <- sample.ras@crs

# extracted cells for which we need covariate values
cells <- extract(sample.ras, sirsam@coords, cellnumbers=TRUE)[,1]

print(paste("Total number of pixels is", 
            as.character(ncell(sample.ras)), ", ", 
            length(cells), 
            "pixels extracted corresponding to shapefile.", sep = " "))

# read covaraites and add to dataframe
intersected.df <- ReadCovariates(covariates.file, cells, sample.crs)

# covariate names in df after intersection
intersected.names <- names(intersected.df)

# regression formula
fmla <- as.formula(paste('Na_ppm_imp', " ~ ", 
                         paste(intersected.names[1:6], collapse= " + ")))
# combine with original df
sirsam@data <- cbind.data.frame(sirsam@data, intersected.df)

# determine bw using croosval
model.covariates.bw <- gwr.sel(
  formula = fmla, data = sirsam, longlat = TRUE
  )

# We choose to specify bw (in degrees)
model.covariates.gauss <- gwr(
  formula = fmla,
  data = sirsam,
  bandwidth = 3.0,
  # coords = sirsam@coords,
  se.fit = TRUE,
  longlat = TRUE,
  hatmatrix = TRUE
)

model.covariates.pred <- gwr(
  formula = fmla,
  data = sirsam,
  bandwidth = 3.0,
  se.fit = TRUE,
  fittedGWRobject = model.covariates.gauss,
  longlat = TRUE,
  fit.points = SpatialPoints(cbind(sirsam$DLong, sirsam$Dlats)),
  # predictions = TRUE,
  cl = cl
)

plot(sirsam$Na_ppm_imp, model.covariates.pred$SDF$pred)
