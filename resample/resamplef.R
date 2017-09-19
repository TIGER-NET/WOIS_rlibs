# Fast nearest neighbour and bilinear resampling
# 
# Author: Stefan Schlaffer
# Date: 21.01.2013
###############################################################################


resamplef <- function(x, y, method = "ngb", filename = "", ...)
{
	if (!method %in% c("ngb","bilinear")) stop("method not supported")
	imethd = switch(method, ngb = 1, bilinear = 2, 2)
	
	naval <- -3.4e+38
	
	if (!canProcessInMemory(y, 3)) {
		if (filename == "") filename <- rasterTmpFile()
	}
	
	# process as a whole
	if (filename == "") {
		
		y <- raster(y)
		rr <- 0.5 + (ymax(x) - yFromRow(y, 1:nrow(y))) / yres(x)
		cc <- (xFromCol(y, 1:ncol(y)) - xmin(x)) / xres(x) + 0.5
		rr1 <- trunc(rr)
		cc1 <- trunc(cc)
		overlap.x <- max(round(xres(y) / xres(x)),1)
		overlap.y <- max(round(yres(y) / yres(x)),1)
		minrow <- max(min(rr1, na.rm=TRUE)-overlap.y, 1)
		maxrow <- min(max(rr1, na.rm=TRUE)+overlap.y, nrow(x))
		mincol <- max(min(cc1, na.rm=TRUE)-overlap.x, 1)
		maxcol <- min(max(cc1, na.rm=TRUE)+overlap.x, ncol(x))
		
		if (maxrow >= minrow & maxcol >= mincol) {				# this is to check if strip is intersecting x
			v <- getValuesBlock(x, minrow, maxrow-minrow+1, 
								   mincol, maxcol-mincol+1,
								   format="matrix")
			
			# prepare target raster
			#u <- matrix(naval, nrow(y), ncol(y))
			v[is.na(v)] <- naval
			
			# image coordinates
			src_iy = na.omit(rr - minrow+1)
			src_ix = na.omit(cc - mincol+1)
			
			res <- .Fortran("resample",
					x = v,
					y = matrix(naval, nrow(y), ncol(y)),
					src_nx = as.integer(ncol(v)),
					src_ny = as.integer(nrow(v)),
					trg_nx = as.integer(ncol(y)),
					trg_ny = as.integer(nrow(y)),
					sig_ix = as.double(src_ix),
					sig_iy = as.double(src_iy),
					NAvalue = naval,
					method = as.integer(imethd),
					DUP = FALSE)
		} else {
			# x does not intersect with strip
			res <- list(y = matrix(naval, nrow(y), ncol(y)))
		}
		
		res$y[res$y == naval] <- NA
		res$y <- t(res$y)
		
		y <- setValues(y, as.vector(res$y))
	
	# process in chunks
	} else {
		
		bs <- blockSize(y)
		pb <- pbCreate(bs$n, label="resample", ...)
		y <- raster(y)
		y <- writeStart(y, filename, ...)
		
		for (i in 1:bs$n) {
			
			# prepare target raster
			u <- matrix(naval, bs$nrows[i], ncol(y))
			
			# read from source raster (let x overlap y min. 1 col/row)
			rr <- 0.5 + (ymax(x) - yFromRow(y, bs$row[i]:sum(bs$nrows[1:i]))) / yres(x)
			cc <- (xFromCol(y, 1:ncol(y)) - xmin(x)) / xres(x) + 0.5
			rr1 <- trunc(rr)
			cc1 <- trunc(cc)
			overlap.x <- max(round(xres(y) / xres(x)),1)
			overlap.y <- max(round(yres(y) / yres(x)),1)
			minrow <- max(min(rr1, na.rm=TRUE)-overlap.y, 1)
			maxrow <- min(max(rr1, na.rm=TRUE)+overlap.y, nrow(x))
			mincol <- max(min(cc1, na.rm=TRUE)-overlap.x, 1)
			maxcol <- min(max(cc1, na.rm=TRUE)+overlap.x, ncol(x))
			
			if (maxrow >= minrow & maxcol >= mincol) {				# this is to check if strip is intersecting x
				v <- getValuesBlock(x, minrow, maxrow-minrow+1, 
									   mincol, maxcol-mincol+1,
									   format="matrix")
				v[is.na(v)] <- naval
		
				# image coordinates
				src_iy = na.omit(rr - minrow+1)
				src_ix = na.omit(cc - mincol+1)
			
				res <- .Fortran("resample",
						x = v,
						y = u,
						src_nx = as.integer(ncol(v)),
						src_ny = as.integer(nrow(v)),
						trg_nx = as.integer(ncol(u)),
						trg_ny = as.integer(nrow(u)),
						sig_ix = as.double(src_ix),
						sig_iy = as.double(src_iy),
						NAvalue = naval,
						method = as.integer(imethd),
						DUP = FALSE)
			} else {
				# x does not intersect with strip
				res <- list(y = u)
			}
			
			res$y[res$y == naval] <- NA
			res$y <- t(res$y)
			
			y <- writeValues(y, as.vector(res$y), bs$row[i])
			pbStep(pb, i)
		}
		
		y <- writeStop(y)
		pbClose(pb)
	}	
	return(y)
}
