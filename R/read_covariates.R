#Parse the input text file containing covaraites
############################################################
ParseText <- function(cov_file) {
  all_cov <- read.table(cov_file)
  len1 = 0
  for(i in 1:length(all_cov$V1)){
    fname <- trimws(as.character(all_cov$V1[i]))
    if(!startsWith(fname , "#")){
      if(len1 == 0)
        fnames <- fname
      else
        fnames <- rbind(fnames,fname)
      len1 <- len1 + 1
    }
  }
  rownames(fnames) <- c()
  return(fnames)
}

#Read input covariates from a text file, intersect with targets or road buffer and return a data frame
############################################################
ReadCovariates <- function(cov_file, to_be_intersected, existing_model = NULL) {
  
  if(!is.null(existing_model)){
    previous_model <- raster(existing_model)
    tmp_w <- values(previous_model)
    tmp_w[is.na(tmp_w)] <- 0
    values(previous_model) <- tmp_w
    previous_model <- values(previous_model) / max(values(previous_model), na.rm = TRUE)
    print(paste("Existing prediction at" , existing_model, "will be used as weight."))
  }
    
  all_cov <- ParseText(cov_file)
  
  print("Reading input covariates.")
  expected_CRS <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  print(paste("Expected CRS is", expected_CRS))
  num_categorical <- 0
  
  for(i in 1:nrow(all_cov)){
    fname <- all_cov[i,]
    print(fname)
    r <- raster(fname)
    num_unique <- length(unique(r[]))
    # if a raster has unique values less than col*row / 10000 is assumed to be categorical
    if(num_unique < round((r@ncols * r@nrows) / 10000)){
      print("categorical")
      num_categorical <- num_categorical + 1
      r[] <- as.factor(r[])
    }
    #CRS check
    if(!compareCRS(expected_CRS, r@crs)){
      print(paste("CRS converted for",fname))
      r@crs <- expected_CRS
    }
    
    #weighting
    if(!is.null(existing_model))
      values(r) <- values(r) * previous_model
    
    #Intersect the covariates with the target locations
    cov_intersected <- data.frame(extract(r, to_be_intersected))
    colnames(cov_intersected) <- r@data@names
    if(i == 1)
      df <- data.frame(cov_intersected)
    else
      df <- cbind.data.frame(df,cov_intersected)
    
    if(r@data@isfactor) 
      df[,i] <- as.factor(df[,i])
    
    # #Convert to categorical
    # for(i in 1:length(cov_stack@layers)) 
    #   if(cov_stack@layers[[i]]@data@isfactor) 
    #     df[,i] <- as.factor(df[,i])
    # #Stack 
    # if(i == 1)
    #   cov_stack <- stack(r)
    # else
    #   cov_stack <- stack(cov_stack, r)
  }
  print(paste(as.character(i) , "covariates were read out of which", as.character(num_categorical), "recognized as categorical covariates."))
  return(df)
}