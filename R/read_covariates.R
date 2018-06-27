#Parse the input text file containing covaraites
############################################################
ParseText <- function(cov_file) {
  all_cov <- read.table(cov_file)
  len1 = 0
  for(i in 1:length(all_cov$V1)){
    fname <- trimws(as.character(all_cov$V1[i]))
    if(!startsWith(fname , "#")){
      if(len1 == 0) {
        fnames <- fname
      } else {
        fnames <- rbind(fnames, fname)
      }
      len1 <- len1 + 1
    }
  }
  rownames(fnames) <- c()
  return(fnames)
}

#Read input covariates from a text file, 
# intersect with targets or road buffer and return a data frame
############################################################
ReadCovariates <- function(cov_file, to_be_intersected, sample.crs) {
  all_cov <- ParseText(cov_file)
  print("Reading input covariates.")
  num_categorical <- 0
  
  for(i in 1:nrow(all_cov)){
    fname <- all_cov[i,]
    print(fname)
    r <- raster(fname)
    num_unique <- length(unique(r[]))
    # if a raster has unique values less than col*row / 10000 is assumed to be categorical
    if(num_unique < round((r@ncols * r@nrows) / 10000)){
      print(paste("categorical covariate found", fname, 'will skip'))
      num_categorical <- num_categorical + 1
      r[] <- as.factor(r[])
      next
    }
    # CRS check and correction
    if(!compareCRS(sample.crs, r@crs)){
      print(paste("CRS converted for",fname))
      r@crs <- sample.crs
    }

    #Intersect the covariates with the target locations
    cov_intersected <- data.frame(extract(r, to_be_intersected))
    colnames(cov_intersected) <- r@data@names
    if(i == 1)
      df <- data.frame(cov_intersected)
    else
      df <- cbind.data.frame(df,cov_intersected)
    
    if(r@data@isfactor) 
      df[,i] <- as.factor(df[,i])
    
  }
  print(paste(as.character(i) , "covariates were read out of which", 
              as.character(num_categorical), 
              "recognized as categorical covariates."))
  return(df)
}
