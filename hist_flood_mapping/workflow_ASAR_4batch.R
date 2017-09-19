#############################################################################
#  
#  Historic Flood mapping workflow for ENVISAT ASAR Wide Swath (WS) data
#
#  Developped in the framework of the ESA-funded project TIGERNET
# 
#  Author: Stefan Schlaffer, Vienna University of Technology
#  Contact: stefan.schlaffer@tuwien.ac.at
#  Date: 18.11.2012
#
#  Depends: R (>= 2.14)
#           R packages raster, rgdal
#		       otsu.R
#
#  Input: Flood image (required)
#         Reference image (optional)
#	  Path and name of the output file (required)
#	  Mask image showing areas prone to flooding (required)
#	  User-defined threshold value
#	  Tile size for splitting image for threshold computation
#
#  Output: Classified GeoTIFF (0: non-flooded, 1: flooded, 250: masked)
#
#############################################################################


library(raster)
library(rgdal)
library(maptools)
library(rgeos)


r1 <- raster(floodimg)
projection(r1) <- CRS("+init=epsg:4326")
error <- try(r1 <- crop(r1,AOI))
if (class(error) == "try-error") {
	cat("no intersect\n")
	next
}
v <- extract(r1, AOI)[[1]]
if (all(is.na(v))) {
	cat("no intersect\n")
	next
}
rm(v)

#r1 <- r1/100
if (refimg != "") useRefimg <- TRUE else useRefimg <- FALSE
if (useRefimg) {
  r2 <- raster(refimg)
  #NAvalue(r2) <- nodatavalue
  r2 <- crop(r2,AOI)
}
fpmask <- raster(floodproneimg)


## collocation
if (useRefimg) {
  if (!compareRaster(r1,r2,stopiffalse = FALSE)) {
    r2 <- resamplef(r2, r1, "bilinear", progress = "window")
  }
}

## output file names
outfile1 <- file.path(workdir, basename(floodimg))
extension(outfile1) <- ""
outfile2 <- file.path(workdir, basename(floodimg))
extension(outfile2) <- ""


## automatic threshold determination
## split-based Otsu
otsuSplitThreshold <- function(img, na.rm = TRUE) {
  ## retrieve image statistics
  imgmean <- cellStats(img, "mean", na.rm=TRUE)
  imgmax <- cellStats(img, "max")
  imgmin <- cellStats(img, "min")
  imggraymean <- (imgmean-imgmin) / abs(imgmax-imgmin) * 255

  ## compute Otsu threshold for each image tile
  nTilesY <- floor(nrow(img) / splitTileSize)
  aspectratio <- round(nrow(img)/ncol(img))
  nTilesX <- floor(ncol(img) / splitTileSize)
  #thr <- vector("numeric", nTilesY * nTilesX)
  tilestats <- as.data.frame(matrix(NA, nr=nTilesY*nTilesX, 6))
  names(tilestats) <- c("i","j","CV","R","nNoData","otsu")

  for (i in 1:nTilesY) {
    cat("row",i,"...\n")
    rowmin <- max(1,(i-1)*splitTileSize+1)
    m <- getValues(img, row = rowmin, nrows = splitTileSize, format = "matrix")

    for (j in 1:nTilesX) {
      colmin <- max(1,(j-1)*splitTileSize/aspectratio+1)
      colmax <- min(ncol(img),(j-1)*splitTileSize/aspectratio+splitTileSize/aspectratio)
      v <- as.vector(m[,colmin:colmax])
      
      if (any(!is.na(v))) {
        v.gray <- (v - imgmin) / (abs(imgmax-imgmin)) * 255
        k <- (i-1)*nTilesX+j
        tilestats[k,1] <- i
        tilestats[k,2] <- j
        tilestats[k,3] <- abs(sd(v.gray, na.rm=TRUE) / mean(v.gray, na.rm=TRUE))
        tilestats[k,4] <- mean(v.gray, na.rm=TRUE)/imggraymean
        tilestats[k,5] <- length(which(is.na(v)))
        if (tilestats[k,5] < splitTileSize^2/100 & tilestats[k,3] >= 0.7 & tilestats[k,4] >= 0.4 & tilestats[k,4] <= 0.9)
          tilestats[k,6] <- otsu(v.gray)
      }
    }
  }

  ## rescale from gray levels
  tilestats$otsu <- tilestats$otsu/255 * (abs(imgmax-imgmin)) + imgmin
  return(tilestats)
}


## classification
cat("Split based threshold computation...\n")
if (useRefimg) {

  r1.lin <- calc(r1, function(x) 10^(x/10))
  r2.lin <- calc(r2, function(x) 10^(x/10))
  changeimg.lin <- r1.lin - r2.lin
  changeimg <- r1 - r2

  tilestats <- try(otsuSplitThreshold(changeimg.lin))
  
  ## use minimum as threshold theta
  if (any(!is.na(tilestats$otsu))) {
    theta <- median(tilestats$otsu, na.rm = TRUE)
    cat("Computed Otsu threshold:", theta, "dB\n")
    r.class.auto <- changeimg.lin < theta
    autoclass <- TRUE
  } else {
    cat("No automatic threshold could be computed.\n")
    autoclass <- FALSE
  }
  r.class.manual <- changeimg < userThreshold
  
} else {

  r1.lin <- calc(r1, function(x) 10^(x/10))
  tilestats <- try(otsuSplitThreshold(r1.lin))
  
  if (class(tilestats) == "try-error") {
    cat("Error occurred during automatic thresholding\n")
    next
  } else {

	  if (any(!is.na(tilestats$otsu))) {
	    
	    ## use minimum as threshold theta
	    theta <- min(tilestats$otsu, na.rm = TRUE)
	    cat("Computed Otsu threshold:", 10*log10(theta), "dB\n")
	    r.class.auto <- r1.lin < theta
	    autoclass <- TRUE
	  } else {
	    cat("No automatic threshold could be computed.\n")
	    autoclass <- FALSE
	  }
	}
  r.class.manual <- r1 < userThreshold
}


## mask non-flood-prone areas
if (!compareRaster(r1,fpmask,stopiffalse = FALSE)) {
  maskfile <- file.path(workdir, basename(floodimg))
  extension(maskfile) <- ""
  maskfile <- paste(maskfile, "_hand_resampled.tif", sep="")
  if (!file.exists(maskfile)) {
    cat("resampling HAND map...\n")
#    fpmask <- resamplef(fpmask, r1, "bilinear", filename = maskfile, datatype="INT2S")
    fpmask <- resamplef(fpmask, r1, "ngb", filename = maskfile, datatype="INT2S")
  } else
    fpmask <- raster(maskfile)
}
cat("masking non-flood-prone areas...\n")
#fpmask <- fpmask < hand.thresh
## combine HAND and user-provided flood-prone areas (UNION):
#fpmask <- Which(fpmask == 1 | rasterize(AOI2, fpmask) == 1)
fpmask[fpmask == 0] <- NA
if (autoclass) r.class.auto <- mask(r.class.auto, fpmask)
r.class.manual <- mask(r.class.manual, fpmask)


## apply MMU
if (mmu > 0) {
  cat("applying MMU...\n")
  if (autoclass) {
    tmp1 <- clump(r.class.auto)
    pxarea <- freq(tmp1)
    r.class.auto[tmp1 %in% pxarea[pxarea[,2] <= mmu,1]] <- 0
  }
  tmp2 <- clump(r.class.manual)
  pxarea <- freq(tmp2)
  r.class.manual[tmp2 %in% pxarea[pxarea[,2] <= mmu,1]] <- 0
}


## output as GeoTIFF
cat("exporting to GeoTIFF...\n")
outfile1 <- paste(outfile1, "_fld_auto.tif", sep="")
outfile2 <- paste(outfile2, "_fld_manual.tif", sep="")
if (autoclass) r.class.auto <- writeRaster(r.class.auto, outfile1, datatype = "INT1U", NAflag = 255)
r.class.manual <- writeRaster(r.class.manual, outfile2, datatype = "INT1U", NAflag = 255)


## output as KMZ
#extension(outfile1) <- "kmz"; extension(outfile2) <- "kmz"
#if (autoclass) KML(r.class.auto, outfile1, maxpixels = ncell(r.class.auto), col = c("grey90","#006699FF"))
#KML(r.class.manual, outfile2, maxpixels = ncell(r.class.manual), col = c("grey90","#006699FF"))
