# Compute Otsu's threshold for a grayscale image
# 
# Author: Stefan Schlaffer, IPF
# Date: 09.03.2011
###############################################################################


otsu <- function(x, na.rm = TRUE)
{
	if (na.rm) x <- na.omit(x)
	if (!is.integer(x)) x <- round(x)
	
	#h <- hist(x, plot = FALSE)
	L <- length(x)
	tt <- table(x)
	vals <- as.integer(names(tt))
	
	wB <- function(p, ct, vals) {
		return(sum(ct[vals < p]) / L)
	}
	wF <- function(p, ct, vals) {
		return(sum(ct[vals >= p]) / L)
	}
	
	vB <- function(p, ct, vals, x) {
		wB(p, ct, vals) * wF(p, ct, vals) * ((mean(x[x < p])-mean(x[x >= p]))^2)
	}
	o <- optimize(vB, c(min(x),max(x)), tt, vals, x, maximum = TRUE)
	return(o$maximum)
}


#.otsu_old <- function(y, m, na.rm = TRUE)
#{
#	if (na.rm) y <- na.omit(y)
#	
#	yvals <- sort(unique(y))
#	
#	L <- length(yvals)
#	per <- as.vector(table(y))/length(y)
#	
#	P <- matrix(0, nrow = L, ncol = L)
#	S <- matrix(0, nrow = L, ncol = L)
#	H <- matrix(0, nrow = L, ncol = L)
#	P[1, ] <- cumsum(per)
#	S[1, ] <- cumsum(per * yvals[1:L])
##	for (u in 2:L) for (v in u:L) {
##		P[u, v] <- P[1, v] - P[1, u - 1]
##		S[u, v] <- S[1, v] - S[1, u - 1]
##	}
#	PS <- .otsu_double_loop(P, S, L)
#	P <- PS[,,1]
#	S <- PS[,,2]
#	rm(PS)
#	
#	H <- S^2/P
#	x <- seq(L)
#	n <- length(x)
#	if (n - 1 < m) 
#		stop("The number of thresholds is larger than the unique values minus 1.")
#	e <- 0
#	h <- m
#	a <- 1:m
#	rule <- c(0, a, L)
#	sigma2 <- sum(sapply(1:(m + 1), function(i) H[rule[i] + 1, 
#								rule[i + 1]]))
#	thresh <- yvals[a]
#	nmmp1 <- n - m + 1
#	mp1 <- m + 1
#	while (a[1] != nmmp1) {
#		if (e < n - h) {
#			h <- 1
#			e <- a[m]
#			j <- 1
#		}
#		else {
#			h <- h + 1
#			e <- a[mp1 - h]
#			j <- 1:h
#		}
#		a[m - h + j] <- e + j
#		if (a[m] != L) {
#			rule <- c(0, a, L)
#			new <- sum(sapply(1:(m + 1), function(i) H[rule[i] + 
#												1, rule[i + 1]]))
#			if (new > sigma2) {
#				sigma2 <- new
#				thresh <- yvals[a]
#			}
#		}
#	}
#	thresh
#}
#
#.otsu2_old <- function(y, m, na.rm = TRUE)
#{
#	if (na.rm) y <- na.omit(y)
#	
#	yvals <- sort(unique(y))
#	
#	L <- length(yvals)
#	per <- as.vector(table(y))/length(y)
#	
#	thresh <- .otsu_f(y, n = length(y),
#			m = m, yvals = yvals,
#			per = per, L = L)
#	return(thresh)
#}
#
#.otsu_double_loop <- function(P, S, L)
#{
#	if (!all(!is.na(P)) | !all(!is.na(S))) {
#		P[is.na(P)] <- -9999
#		S[is.na(S)] <- -9999
#	}
#	
#	PS <- array(as.double(0), dim = c(L, L, 2))
#	PS[,,1] <- P
#	PS[,,2] <- S
#	
#	PS <- .Fortran("otsu_double_loop",
#			PS = PS, L = as.integer(L),
#			DUP = FALSE)
#	return(PS$PS)
#}
#
#
#.otsu_f <- function(y, n, m, yvals, per, L)
#{
#	thresh <- rep(0, m)
#	
#	out <- .Fortran("otsu",
#			y = as.double(y),
#			n = as.integer(n),
#			m = as.integer(m),
#			yvals = as.double(yvals),
#			per = as.double(per),
#			L = as.integer(L),
#			thresh = as.double(thresh))
#	
#	return(out$thresh)
#}
