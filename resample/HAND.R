# Height Above Nearest Drainage (HAND) index as 
# described by Renno et al. (2008)
# 
# Author: Stefan Schlaffer
# Date: 13.07.2012
# TODO: make safe for large files
###############################################################################


HAND <- function(dem, flowdir, drainage, filename = "", flowdir.type = "ArcGIS", ...)
{
	if (!canProcessInMemory(dem,1))
		stop("Data too large to process in memory. Try a subset.")
	if (!is.integer(head(flowdir)) | !is.integer(head(drainage)))
		stop("flowdir and drainage must be INTEGER data types")
	
	compareRaster(dem, flowdir, drainage)
	
	demdatatype <- ifelse(is.integer(head(dem)), "INT", "FLT")
	
	# label drainage network cells
	cat("preparing drainage network...\n")
	draincells <- Which(drainage != 0 & !is.na(drainage), cells = TRUE)
	drainage[draincells] <- 1:length(draincells)
	if (!flowdir.type %in% c("ArcGIS","GRASS"))
		stop("flowdir.type must be one of c('ArcGIS','GRASS')")
	grass.flowdir <- flowdir.type == "GRASS"
	
	# pass matrices to Fortran function
	cat("labelling cells...\n")
	nr <- nrow(flowdir)
	nc <- ncol(flowdir)
	flowdir.m <- as.matrix(flowdir)
	drainage.m <- as.matrix(drainage)
	dem.m <- as.matrix(dem)
	flowdir.m[is.na(flowdir.m)] <- as.integer(-32768)
	drainage.m[is.na(drainage.m)] <- as.integer(-32768)
	if (demdatatype == "INT") 
		dem.m[is.na(dem.m)] <- as.integer(-32768)
	else if (demdatatype == "FLT")
		dem.m[is.na(dem.m)] <- as.numeric(-32768)
	#ht.nr <- nrow(ht)
	
	if (demdatatype == "INT") {
		fout <- .Fortran("hand",
				dir = flowdir.m,
				drain = drainage.m,
				dem = dem.m,
				nr = as.integer(nr),
				nc = as.integer(nc),
				#ht = ht,
				#ht_nr = as.integer(ht.nr),
				res = matrix(as.integer(0), nr, nc),
				nodataval = as.integer(-32768),
				dir_grass = as.integer(grass.flowdir),
				DUP = FALSE)
	} else if (demdatatype == "FLT") {
		fout <- .Fortran("hand_flt",
				dir = flowdir.m,
				drain = drainage.m,
				dem = dem.m,
				nr = as.integer(nr),
				nc = as.integer(nc),
				#lab = as.integer(ht[,1]),
				#ht = as.numeric(ht[,2]),
				#ht_nr = as.integer(ht.nr),
				res = matrix(as.numeric(0), nr, nc),
				nodataval = as.numeric(-32768),
				dir_grass = as.integer(grass.flowdir),
				DUP = FALSE)
	}
	
	fout$res[fout$res == -32768] <- NA
	fout$res[fout$res < 0] <- 0
	
	# output
	res <- raster(dem)
	if (!canProcessInMemory(dem, 4)) {
		if (filename == "") filename <- rasterTmpFile()
		bs <- blockSize(res)
		res <- writeStart(res, filename, ...)
		fout$res <- t(fout$res)
		for (i in 1:bs$n) {
			res <- writeValues(res, fout$res[((bs$row[i]-1)*ncol(res)+1):(ncol(res)*(sum(bs$nrows[1:i])))], bs$row[i])
		}
		res <- writeStop(res)
	} else { 
		res[] <- fout$res
		if (filename != "") {
			res <- writeRaster(res, filename, ...)
		}
	}
	
	return(res)
}
