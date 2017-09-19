# Apply a Minimum Mapping Unit using the spatstat package
# 
# Author: IPF\ss
###############################################################################


applyMMU <- function(x, size, filename = "", ...)
{
	if (class(x) != "RasterLayer") {
		stop("x has to be of class RasterLayer")
	}
	
	require(spatstat)
	
	# always do chunkwise processing
#	if (!canProcessInMemory(x, 2) & filename == "")
	filename <- rasterTmpFile()
	
	if (filename == "") {
		v <- getValues(x, format="matrix")
		if (any(is.finite(v), na.rm=TRUE)) {
		
			if (any(v > 0)) {
				u <- v
				u[u == 0] <- NA
				con <- connected(im(u), background=0)
				con <- as.matrix(con)
				
				cnt <- table(con)
				cnt <- as.data.frame(cnt)
				v[which(con %in% cnt[cnt$Freq <= size,1])] <- 0
#				ind.rm <- which(as.matrix(cnt) <= size)
#				if (length(ind.rm) > 0) {
#					for (k in 1:length(ind.rm)) {
#						v[which(con == ind.rm[k])] <- 0
#					}
#				}
			}
		}
		
		out <- raster(x)
		out <- setValues(out, t(v))
		
	} else {
		bs <- blockSize(x)
		pb <- pbCreate(bs$n, ...)
		out <- raster(x)
		out <- writeStart(out, filename, ...)
		for (i in 1:bs$n) {
			rowmin <- max(bs$row[i]-size, 1)
			nrows <- min(bs$nrows[i]+2*size, nrow(x)-rowmin)
			v <- getValues(x, rowmin, nrows, format="matrix")
			
			if (any(is.finite(v))) {
			
				if (any(v > 0, na.rm=TRUE)) {
					u <- v
					ina <- which(is.na(u))
					u[u == 0] <- NA
					con <- connected(im(u), background=NA)
					con <- as.matrix(con)
					
					cnt <- table(con)
					cnt <- as.data.frame(cnt)
					v[which(con %in% cnt[cnt$Freq <= size,1])] <- 0
					v[ina] <- NA
#					ind.rm <- which(as.matrix(cnt) <= size)
#					if (length(ind.rm) > 0) {
#						for (k in 1:length(ind.rm)) {
#							v[which(con == ind.rm[k])] <- 0
#						}
#					}
				}
			}
			if (i != 1 & i != bs$n) {
				v <- v[(size+1):(size+bs$nrows[i]),]
			} else if (i == 1) {
				v <- v[1:bs$nrows[i],]
			} else {
				v <- v[(size):nrow(v),]
			}
			out <- writeValues(out, t(v), bs$row[i])
			pbStep(pb, i)
		}
		out <- writeStop(out)
		pbClose(pb, i)
	}
	
	return(out)
}

