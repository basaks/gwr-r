rm(list=ls())
library(sp)
library(rgdal)
library(GWmodel)
library(doParallel)
library(raster)

source('R/read_covariates.R')
sirsam <- readOGR(dsn = "/g/data1/ge3/john/jobs/murray_B_depth/Targets_2017/Targets_MurrayB_confid3_resampled_25p.shp")
covariates.file <- 'list.txt'
small.covariates.file <- 'list.small.txt'

LONGLAT <- FALSE
ADAPTIVE <- TRUE
KERNEL <- 'gaussian'
APPROACH <- 'AIC'

# pick a sample raster for gathering the indices of the targets in rasters
covariates.list <- ParseText(covariates.file)
sample.ras <- raster(covariates.list[1])
sample.crs <- sample.ras@crs


small.covariates.list <- ParseText(small.covariates.file)
small.sample.ras <- raster(small.covariates.list[1])
small.sample.crs <- small.sample.ras@crs


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
fmla <- as.formula(paste('depth', " ~ ",
                         paste(intersected.names, collapse= " + ")))

# combine with original df
sirsam@data <- cbind.data.frame(sirsam@data, intersected.df)

# determine bw using croosval
model.covariates.bw <- bw.gwr(
  formula = fmla, data = sirsam,
  longlat = LONGLAT, kernel = KERNEL,
  approach = APPROACH, adaptive = ADAPTIVE
)

# We choose to specify bw (in degrees)
# model.covariates.gauss <- gwr.basic(
#   formula = fmla,
#   data = sirsam,
#   bw = model.covariates.bw,
#   # coords = sirsam@coords,
#   longlat = LONGLAT, kernel = KERNEL,
#   adaptive = ADAPTIVE
# )

model.covariates.pred <- gwr.predict(
  formula = fmla,
  data = sirsam,
  bw = model.covariates.bw,
  predictdata = SpatialPointsDataFrame(
    coords = coordinates(sirsam),
    data = sirsam@data,
    proj4string = sample.crs),
  longlat = LONGLAT, kernel = KERNEL,
  adaptive = ADAPTIVE
)

# print measured_vs_prediction output
pdf('measured_vs_prediction.pdf')
plot(sirsam$depth, model.covariates.pred$SDF$prediction,
     main="Prediction at Mesurement Points",
     xlab="Measured", ylab="Prediction")
 dev.off()
# plot(sirsam$Na_ppm_imp, sqrt(model.covariates.pred$SDF$prediction_var))

predictors <- stack()

for (i in 1:nrow(small.covariates.list)) {
  predictors <- stack(predictors, small.covariates.list[i])
}

row.max <- dim(predictors)[1]
col.max <- dim(predictors)[2]
col.num <- 0
all.coords <- coordinates(predictors)


# precalculate sysmetric distance matrix between data points
dm.calib <- gw.dist(dp.locat = coordinates(sirsam))

predfunc <- function(i){
  library(GWmodel)
  library(raster)
  print(paste('predicting for loop: ', i))
  last.row <- i + 1
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
    longlat = LONGLAT,
    kernel = KERNEL,
    adaptive = ADAPTIVE,
    dMat2 = dm.calib
  )
  return(chunk.pred$SDF$prediction)
}

# run the parallel computation
cl <- makeCluster(detectCores())
registerDoParallel(cl)
finalRasterPred <- foreach(i=1:row.max) %dopar% predfunc(i)
stopCluster(cl)

s2 <- writeStart(small.sample.ras, filename='mb_prediction.tif', format='GTiff',
                 overwrite=TRUE, dataType='FLT4S')

for (i in seq(from=1, to=row.max, by=1)) {
  v <- unlist(finalRasterPred[i])
  s2 <- writeValues(s2, v, i)
}
s2 <- writeStop(s2)
