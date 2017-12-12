rm(list=ls())
library(sp)
library(rgdal)
library(GWmodel)
library(parallel)
library(raster)


source('R/read_covariates.R')
sirsam <- readOGR(dsn = "data/geochem_sites.shp")
covariates.file <- 'data/sirsam_covariates_Na.txt'

LONGLAT <- TRUE
ADAPTIVE <- TRUE
KERNEL <- 'gaussian'
APPROACH <- 'AIC'


cl <- makeCluster(detectCores() - 1)


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
model.covariates.bw <- bw.gwr(
  formula = fmla, data = sirsam, 
  longlat = LONGLAT, kernel = KERNEL,
  approach = APPROACH, adaptive = ADAPTIVE
)

# We choose to specify bw (in degrees)
model.covariates.gauss <- gwr.basic(
  formula = fmla,
  data = sirsam,
  bw = model.covariates.bw,
  # coords = sirsam@coords,
  longlat = LONGLAT, kernel = KERNEL,
  approach = APPROACH, adaptive = ADAPTIVE
)


model.covariates.pred <- gwr.predict(
  formula = fmla,
  data = sirsam,
  bw = model.covariates.bw,
  predictdata = SpatialPointsDataFrame(
    coords = coordinates(sirsam),
    data = sirsam@data, 
    proj4string = sample.crs),
  longlat = LONGLAT, kernel = KERNEL,
  approach = APPROACH, adaptive = ADAPTIVE
)

plot(sirsam$Na_ppm_imp, model.covariates.pred$SDF$prediction)
plot(sirsam$Na_ppm_imp, sqrt(model.covariates.pred$SDF$prediction_var))

predictors <- stack()

for (i in 1:nrow(covariates.list)) {
  predictors <- stack(predictors, covariates.list[i])
}

predfun <- function(model.covariates.gauss, spatial.df){
  v <- gwr.predict(
    formula = fmla,
    data = sirsam,
    bw = model.covariates.bw,
    predictdata = spatial.df,
    longlat = TRUE,
    kernel = 'gaussian',
    adaptive = TRUE
  )
  return(cbind(p=as.vector(v$SDF$prediction), 
               var=as.vector(v$SDF$prediction_var)))
}

row.block.size <- 2
row.max <- dim(predictors)[1]
col.max <- dim(predictors)[2]
col.num <- 0
all.coords <- coordinates(predictors)


s2 <- writeStart(sample.ras, filename='test2.tif', format='GTiff', 
                 overwrite=TRUE)


for (i in seq(from=1, to=row.max, by=row.block.size)) {
  last.row <- i + row.block.size
  if (last.row > row.max) {
    last.row <- row.max +1
  }
  
  start.i <- (i-1)*col.max + 1
  end.i <- (last.row-1)*col.max
  
  coords <- all.coords[start.i: end.i, ]
  last.p <- last.row-1
  df <- data.frame(predictors[i: last.p, ])
  
  predictdata <- SpatialPointsDataFrame(
    coords = coords,
    data = df,
    proj4string = sample.crs
    )

  chunk.pred <- gwr.predict(
    formula = fmla,
    data = sirsam,
    bw = model.covariates.bw,
    predictdata = predictdata,
    longlat = TRUE,
    kernel = 'gaussian',
    adaptive = TRUE
  )
  
  s2 <- writeValues(s2, chunk.pred$SDF$prediction, i)
  
  col.num <- col.num + 1
  coords <- 0
  }
}
s2 <- writeStop(s2)
