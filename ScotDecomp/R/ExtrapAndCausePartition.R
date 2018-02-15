
# Author: tim
###############################################################################

# a script to extrapolate mx from age 85+ to age 110+
# header to determine working directory
me <- system("whoami",intern=TRUE)
if (me == "tim"){
	setwd("/home/tim/git/ScotDecomp/ScotDecomp")
}
if (me == "mpidr_d\\seaman"){
	setwd("U:/Conferences/PAA/2018 Denver/within and between/ScotDecomp")
	
}

# -----------------------------------------------
source("R/Functions.R")
library(devtools)

# this can be modified if necessary.
#source("/home/tim/git/DistributionTTD/DistributionTTD/R/Functions.R")
library(ungroup)
#install_github("mpascariu/ungroup", dependencies = TRUE)
# install.packages("MortalitySmooth")
library(MortalitySmooth)
library(reshape2)
library(HMDHFDplus)
library(data.table)
#install_github("timriffe/DecompHoriuchi/DecompHoriuchi", dependencies = TRUE)
library(DecompHoriuchi)
library(RColorBrewer)

graduate <- function(chunk){
    chunk      <- as.data.frame(chunk)
	chunkold   <- data.frame(age=100:110)
	chunkold   <- merge(chunkold, chunk, all.x = TRUE)
	
	imputethingy <- function(x){
		if (all(is.na(x))){
			return(x)
		}
		unique(na.omit(x))
	}

	chunkold$Fac      <- imputethingy(chunk$Fac)
	chunkold$quintile <- imputethingy(chunk$quintile)
	chunkold$year     <- imputethingy(chunk$year)
	chunkold$sex      <- imputethingy(chunk$sex)
	
	chunk      <- rbind(chunk[1:100, ], chunkold[, colnames(chunk)])
	
	chunk2     <- chunk[!is.na(chunk$N),]
	chunk2$age <- age2int(chunk2$age)
	x          <- chunk2$age
	y          <- chunk2$N
	

	nlast      <- 111 - max(x)
	TheNewN    <- pclm(x, y, nlast, out.step = 1,show=FALSE)
	chunk$N2   <- TheNewN$fitted
	chunk[order(chunk$age), ]
}

extrap <- function(DX, NX, fit.age.i = 81, extrap.age.i = 81){
	extrap.ages.i <- rep(extrap.age.i,ncol(DX))
	suppressWarnings(ltper_mx_v5(DX,NX,extrap.ages.i=extrap.ages.i,fit.from.i = fit.age.i ))
}

SmoothFracs <- function(chunk){
	chunk      <- chunk[order(chunk$age), ]
	Dxc        <- as.matrix(chunk[,c("agg_other", "agg_external", "agg_amn", "agg_resp","agg_circular")])
	FRACS      <- Dxc * 0
	colnames(FRACS) <- paste0(colnames(Dxc)[1:5],"_frac_sm")
	m0         <- Dxc[1,1:5]/chunk$N3[1]
	FRACS[1, ] <- m0
	x          <- 1:110
	EXP        <- chunk$N3 * chunk$Fac
	offset     <- log(EXP[-1])
	for (i in 1:5){
		y                  <- Dxc[-1,i]
		FRACS[2:111, i]    <- exp(Mort1Dsmooth(x, y, offset, method=3, lambda=1)$logmortality)
	}

	#y          <- chunk$agg_allc[-1] * chunk$Fac[-1]
	y          <- chunk$mx * chunk$N3
	mxsm       <- exp(Mort1Dsmooth(x, y[-1], 
					offset = log(chunk$N3)[-1], method = 3, lambda = .1)$logmortality)
	mxsm       <- c(chunk$mx[1], mxsm)
	chunk$mxsm <- mxsm
#	plot(0:110, chunk$mx,type='l',log='y')
#	lines(0:110,mxsm)
#	matplot(FRACS,type='l',log='y')
#	matplot(as.matrix(chunk[,1:5])/(chunk$N3*3),type='l',log='y',add=TRUE)
	FRACS      <- FRACS / rowSums(FRACS)
	FRACS      <- as.data.frame(FRACS)
	chunk      <- cbind(chunk, FRACS)
	
	FRACS      <- as.data.frame(Dxc / rowSums(Dxc))
	colnames(FRACS) <-paste0(colnames(FRACS),"_frac")
	chunk      <- cbind(chunk, FRACS)
	chunk
}
# ---------------------------------
SCO             <- read.csv("Data/Side/deaths_Scotland.csv")
# remove some chaff columns
SCO$population  <- NULL
SCO$popcheck    <- NULL
SCO$deathcheck  <- NULL

# rebase pops and denoms, no real reason
SCO$Fac                   <- 3
SCO$Fac[SCO$year == 1991] <- 2
SCO$N                     <- SCO$N / SCO$Fac
SCO$agg_allc              <- SCO$agg_allc / SCO$Fac
# chop off last year
SCO                       <- SCO[SCO$age <= 110, ]



# Graduate, split, apply, combine
SCOL            <- split(SCO, list(SCO$year,SCO$sex,SCO$quintile))
SCOL            <- lapply(SCOL, graduate)
SCO             <- do.call(rbind,SCOL)
SCO             <- SCO[order(SCO$year, SCO$sex, SCO$quintile, SCO$age), ]
# single ages altered. Let's only keep the graduated 
# data points that we actually needed
# and leave single ages alone
SCO$N3          <- SCO$N
ind             <- is.na(SCO$N) | SCO$age_cat != ""
SCO$N3[ind]     <- SCO$N2[ind]

# no NAs, only 0s
SCO$agg_allc[is.na(SCO$agg_allc)] <- 0

# now to extrapolate to 110, create data arrays
DXm             <- acast(SCO[SCO$sex == 1, ], age~year~quintile, value.var = "agg_allc")
DXf             <- acast(SCO[SCO$sex == 2, ], age~year~quintile, value.var = "agg_allc")
NXm             <- acast(SCO[SCO$sex == 1, ], age~year~quintile, value.var = "N3")
NXf             <- acast(SCO[SCO$sex == 2, ], age~year~quintile, value.var = "N3")

DXm[is.na(DXm)] <- 0
DXf[is.na(DXf)] <- 0

# this is where different versions of results will go
# the first number is the lowest age used for information
# the second number is the age from which we extrapolate.
# I chose 75.85 based on visual inspection.

Container.mx.m.70.80 <- array(dim=c(111,dim(DXm)[2:3]),dimnames=list(0:110,seq(1981,2011,10),1:5))
Container.mx.f.70.80 <- array(dim=c(111,dim(DXf)[2:3]),dimnames=list(0:110,seq(1981,2011,10),1:5))
Container.mx.m.75.80 <- array(dim=c(111,dim(DXm)[2:3]),dimnames=list(0:110,seq(1981,2011,10),1:5))
Container.mx.f.75.80 <- array(dim=c(111,dim(DXf)[2:3]),dimnames=list(0:110,seq(1981,2011,10),1:5))
Container.mx.m.80.80 <- array(dim=c(111,dim(DXm)[2:3]),dimnames=list(0:110,seq(1981,2011,10),1:5))
Container.mx.f.80.80 <- array(dim=c(111,dim(DXf)[2:3]),dimnames=list(0:110,seq(1981,2011,10),1:5))

Container.mx.m.70.85 <- array(dim=c(111,dim(DXm)[2:3]),dimnames=list(0:110,seq(1981,2011,10),1:5))
Container.mx.f.70.85 <- array(dim=c(111,dim(DXf)[2:3]),dimnames=list(0:110,seq(1981,2011,10),1:5))
Container.mx.m.75.85 <- array(dim=c(111,dim(DXm)[2:3]),dimnames=list(0:110,seq(1981,2011,10),1:5))
Container.mx.f.75.85 <- array(dim=c(111,dim(DXf)[2:3]),dimnames=list(0:110,seq(1981,2011,10),1:5))
Container.mx.m.80.85 <- array(dim=c(111,dim(DXm)[2:3]),dimnames=list(0:110,seq(1981,2011,10),1:5))
Container.mx.f.80.85 <- array(dim=c(111,dim(DXf)[2:3]),dimnames=list(0:110,seq(1981,2011,10),1:5))


# just ignore warnings and printouts
for (i in 1:dim(DXm)[3]){
	DXmi <- DXm[,,i]
	DXfi <- DXf[,,i]
	NXmi <- NXm[,,i]
	NXfi <- NXf[,,i]
	DXmi[96:111,] <- NA
	DXfi[96:111,] <- NA
	Container.mx.m.70.80[,,i] <- extrap(DXmi,NXmi,fit.age.i =71,extrap.age.i=81)
	Container.mx.f.70.80[,,i] <- extrap(DXfi,NXfi,fit.age.i =71,extrap.age.i=81)
	Container.mx.m.75.80[,,i] <- extrap(DXmi,NXmi,fit.age.i =76,extrap.age.i=81)
	Container.mx.f.75.80[,,i] <- extrap(DXfi,NXfi,fit.age.i =76,extrap.age.i=81)
	Container.mx.m.80.80[,,i] <- extrap(DXmi,NXmi,fit.age.i =81,extrap.age.i=81)
	Container.mx.f.80.80[,,i] <- extrap(DXfi,NXfi,fit.age.i =81,extrap.age.i=81)
	Container.mx.m.70.85[,,i] <- extrap(DXmi,NXmi,fit.age.i =71,extrap.age.i=86)
	Container.mx.f.70.85[,,i] <- extrap(DXfi,NXfi,fit.age.i =71,extrap.age.i=86)
	Container.mx.m.75.85[,,i] <- extrap(DXmi,NXmi,fit.age.i =76,extrap.age.i=86)
	Container.mx.f.75.85[,,i] <- extrap(DXfi,NXfi,fit.age.i =76,extrap.age.i=86)
	Container.mx.m.80.85[,,i] <- extrap(DXmi,NXmi,fit.age.i =81,extrap.age.i=86)
	Container.mx.f.80.85[,,i] <- extrap(DXfi,NXfi,fit.age.i =81,extrap.age.i=86)
}

mx2edagHMD <- compiler::cmpfun(function(mx, sex = "m"){
			mx                  <- Mna0(as.numeric(mx))
			
			# mean proportion of interval passed at death
			ax                  <- mx * 0 + .5                      # ax = .5, pg 38 MPv5
			
			ax[1]               <- AKm02a0(mx[1], sex)
			
			qx                  <- mx / (1 + (1 - ax) * mx)          # Eq 60 MPv5 (identity)
# ---------------------------------------------------------------------------------
# set open age qx to 1
			i.openage           <- length(mx) # removed argument OPENAGE
			qx[i.openage]       <- 1
			ax[i.openage]       <- 1 / mx[i.openage-1]                   
# ---------------------------------------------------------------------------------
# define remaining lifetable columns:
			px                  <- 1 - qx                                                                                 # Eq 64 MPv5
			px[is.nan(px)]      <- 0 # skips BEL NAs, as these are distinct from NaNs
# lx needs to be done columnwise over px, argument 2 refers to the margin.
			lx                  <- c(1, cumprod(px[1:(i.openage-1)]))
			# NA should only be possible if there was a death with no Exp below age 80- impossible, but just to be sure
			# lx[is.na(lx)]   <- 0 # removed for BEL testing        
			dx                  <- lx * qx                                                                                # Eq 66 MPv5
			
			Lx                  <- lx - (1 - ax) * dx                                                         # Eq 67 MPv5
			Lx[i.openage]     <- lx[i.openage] * ax[i.openage]
# we need to do operations on Lx, but taking its NAs to mean 0
# Lx[is.na(Lx)]    <- 0 # removed for BEL testing
# Tx needs to be done columnwise over Lx, argument 2 refers to the column margin.
			Tx                      <- c(rev(cumsum(rev(Lx[1:(i.openage-1)]))),0) + Lx[i.openage] # Eq 68 MPv5
			ex                      <- Tx / lx 
			# ad hoc return
			dx                      <- dx / sum(dx)
			sum(ex * dx)
		})
mx2edagHMD_vec <- function(mxc,sex = "m"){
	dim(mxc) <- c(111, 5)
	mx       <- rowSums(mxc)
	mx2edagHMD(mx, sex = sex)
}
# chunk here just age and sex, not quintile, we select out quintiles

do.decomp <- function(chunk, version = "sm", quintile1 = 1, quintile2 = 5){
	
	base_names <- c("mx_other", "mx_external", "mx_amn", "mx_resp", "mx_circ")
	if (version == "sm"){
		base_names <- paste0(base_names, "_sm")
	}
	chunk  <- chunk[order(chunk$quintile, chunk$age), ]
	q1     <- as.matrix(chunk[chunk$quintile == quintile1, base_names])
	q2     <- as.matrix(chunk[chunk$quintile == quintile2, base_names])
	
	
	rates1 <- c(q1)
	rates2 <- c(q2)
	
	sex    <- unique(chunk$sex)
	sex    <- ifelse(sex == 1, "m", "f")
	
	# decompose
	dec <- DecompContinuousOrig(
			mx2edagHMD_vec, 
			rates1 = rates1, 
			rates2 = rates2, 
			sex = sex, 
			N = 20)
	# housekeeping	  
	dim(dec)           <- c(111,5)	  
	dimnames(dec)      <- list(0:110, base_names)
	
	year               <- unique(chunk$year)
	attr(dec, "sex")   <- sex
	attr(dec, "year")  <- year
	dec
}

plotdec <- function(X, col = brewer.pal(5, "Dark2"), ylim){
	sex      <- attributes(X)$sex
	year     <- attributes(X)$year
	colns    <- dimnames(X)[[2]]
	smoothed <- grepl("sm",colns[1])
	causes   <- unlist(lapply(strsplit(colns,split="_"),"[[",2))
	
	neg      <- pos <- t(X)
	neg[neg > 0] <- 0
	pos[pos < 0] <- 0
	
	if (missing(ylim)){
		ylim <- c(min(colSums(neg)), max(colSums(pos)))
	}
	
	xgrid <-  seq(0, 100, by = 10)
	barplot(pos, width = rep(1,111), space = 0, xlim = c(0, 100), 
			ylim = ylim, axes = FALSE, border = NA, 
			axisnames = FALSE, legend.text = causes,
			main = paste(year,sex,ifelse(smoothed,"smoothed","")),col=col)
	par(new=TRUE)
	barplot(neg, width = rep(1,111), space = 0, xlim = c(0, 100), ylim=ylim,
			add = TRUE, axisnames = FALSE, axes = FALSE, border = NA,col=col)
	axis(1, pos = ylim[1] * 1.05)
	segments( xgrid, ylim[1] * 1.05, xgrid, .1, col = "white", lwd = .7)
	axis(2, las = 1, at = pretty(ylim)
	NULL
}


# all results generated, but we only merge in 75.85
mxmlong     <- melt(Container.mx.m.75.85, value.name="mx",varnames=c("age","year","quintile"))
mxflong     <- melt(Container.mx.f.75.85, value.name="mx",varnames=c("age","year","quintile"))
mxmlong$sex <- 1
mxflong$sex <- 2
# combine into single object
SCOlong     <- rbind(mxmlong, mxflong)
SCO         <- merge(SCO, SCOlong)
SCO         <- SCO[order(SCO$year, SCO$sex, SCO$quintile, SCO$age), ]


# ----------------------------------------------
# now for cause fractions
# derive COD fractions from smoothed counts, split, apply, combine, merge
SCOL               <- split(SCO, list(SCO$year, SCO$sex, SCO$quintile))
SCOL               <- lapply(SCOL, SmoothFracs)
SCO                <- do.call(rbind, SCOL)


# get Mxc
SCO$mx_other       <- SCO$agg_other_frac_sm * SCO$mx
SCO$mx_external    <- SCO$agg_external_frac_sm * SCO$mx
SCO$mx_amn         <- SCO$agg_amn_frac_sm * SCO$mx
SCO$mx_resp        <- SCO$agg_resp_frac_sm * SCO$mx
SCO$mx_circ        <- SCO$agg_circular_frac_sm * SCO$mx

# create the same from mxsm 
SCO$mx_other_sm    <- SCO$agg_other_frac_sm * SCO$mxsm
SCO$mx_external_sm <- SCO$agg_external_frac_sm * SCO$mxsm
SCO$mx_amn_sm      <- SCO$agg_amn_frac_sm * SCO$mxsm
SCO$mx_resp_sm     <- SCO$agg_resp_frac_sm * SCO$mxsm
SCO$mx_circ_sm     <- SCO$agg_circular_frac_sm * SCO$mxsm

# visual diagnostic
SCOL               <-  split(SCO, list(SCO$year, SCO$sex, SCO$quintile))

# start here, make sure all mxc are normal looking.
# plot raw Mxc, semismooth Mxc, and smooth Mxc
pdf("Figures/Diagnostics/AllCauseSmoothing.pdf")
lapply(SCOL, function(chunk){
			sex   <- unique(chunk$sex)
			sex   <- ifelse(sex == 1, "male","female")
			yr    <- unique(chunk$year)
			quint <- unique(chunk$quintile)
		plot(0:110, chunk$mx, pch = 16, cex = .7, log = 'y', main = paste(sex, yr, "q=",quint),
				xlab = "Age",ylab = "log rate",ylim=c(1e-7,1))
		lines(0:110, chunk$mxsm)
		})
dev.off()


pdf("Figures/Diagnostics/MxcCompare.pdf")
lapply(SCOL, function(chunk){
			sex   <- unique(chunk$sex)
			sex   <- ifelse(sex == 1, "male","female")
			yr    <- unique(chunk$year)
			quint <- unique(chunk$quintile)
			
			cols  <- c("agg_other", "agg_external", "agg_amn", "agg_resp", "agg_circular")
			Mxc   <- as.matrix(chunk[,cols]) / (chunk$N3*chunk$Fac)
			
			cols  <- c("mx_other", "mx_external", "mx_amn", "mx_resp", "mx_circ")    
			mxc   <- as.matrix(chunk[,cols])
			
			cols  <- c("mx_other_sm", "mx_external_sm", "mx_amn_sm", "mx_resp_sm", "mx_circ_sm")   
			mxcsm <- as.matrix(chunk[,cols])
			
			matplot(0:110, Mxc, pch = 16, cex = .7, log = 'y', main = paste(sex, yr, "q=",quint),
					xlab = "Age", ylab = "log rate", ylim=c(1e-7,1))
			matplot(0:110, mxc, type = 'l', lwd = .8, add = TRUE)
			matplot(0:110, mxcsm, type = 'l', lwd = 1.2, add = TRUE)
		})
dev.off()

pdf("Figures/Diagnostics/Mxcsm.pdf")
lapply(SCOL, function(chunk){
			sex   <- unique(chunk$sex)
			sex   <- ifelse(sex == 1, "male","female")
			yr    <- unique(chunk$year)
			quint <- unique(chunk$quintile)

			cols  <- c("mx_other_sm", "mx_external_sm", "mx_amn_sm", "mx_resp_sm", "mx_circ_sm")   
			mxcsm <- as.matrix(chunk[,cols])
			
			matplot(0:110, mxcsm,type = 'l',  lwd = 1.2, log = 'y', main = paste(sex, yr, "q=",quint),
					xlab = "Age", ylab = "log rate", ylim=c(1e-7,1))
		})
dev.off()

# Now decompose and plot
SCO <- as.data.frame(SCO)
SCOL   <- split(SCO, list(SCO$year, SCO$sex))
DECsmL <- lapply(SCOL, do.decomp)
DECL   <- lapply(SCOL, do.decomp, version = "")

cols <- brewer.pal(5,"Dark2")
pdf("Figures/Diagnostics/DecompSm.pdf")
lapply(DECsmL, plotdec, col = cols, ylim = c(-0.1,.35))
dev.off()

pdf("Figures/Diagnostics/Decomp.pdf")
lapply(DECL, plotdec, col = cols, ylim = c(-0.1,.35))
dev.off()

SCOL   <- split(SCO, list(SCO$year, SCO$sex, SCO$quintile))
pdf("Figures/Diagnostics/MxDiffs.pdf")
lapply(SCOL, function(X){
			mx     <- X$mx
			mxsm   <- X$mxsm
			dmx    <- diff(log(mx))[-1]
			dmxsm  <- diff(log(mxsm))[-1]
			
			ylim   <- range(pretty(c(dmx,dmxsm)))
			plot(1:109, dmx, type = 'l', ylim = ylim)
			lines(1:109, dmxsm)
		})
dev.off()
SCOL   <- split(SCO, list(SCO$year, SCO$sex))
# look at quintile difference in log rate schedules all cause and circ,
# raw and smoothed.
pdf("Figures/Diagnostics/MxDiffs1v5.pdf")
lapply(SCOL, function(X){
			
			mx1        <- X$mx[X$quintile == 1]
			mx5        <- X$mx[X$quintile == 5]
			mx1sm      <- X$mxsm[X$quintile == 1]
			mx5sm      <- X$mxsm[X$quintile == 5]
			
			mx1c        <- X$mx_circ[X$quintile == 1]
			mx5c        <- X$mx_circ[X$quintile == 5]
			mx1smc      <- X$mx_circ_sm[X$quintile == 1]
			mx5smc      <- X$mx_circ_sm[X$quintile == 5]
			
			mx1a        <- X$mx_amn[X$quintile == 1]
			mx5a        <- X$mx_amn[X$quintile == 5]
			mx1sma      <- X$mx_amn_sm[X$quintile == 1]
			mx5sma      <- X$mx_amn_sm[X$quintile == 5]
			
			
			diffs       <- log(mx1) - log(mx5)
			diffssm     <- log(mx1sm) - log(mx5sm)
			
			diffsc      <- log(mx1c) - log(mx5c)
			diffssmc    <- log(mx1smc) - log(mx5smc)
			
			diffsa      <- log(mx1a) - log(mx5a)
			diffssma    <- log(mx1sma) - log(mx5sma)
			
			ylim       <- range(pretty(c(diffs,diffsc)))
			plot(0:110, diffs, type = 'l', ylim = ylim)
			lines(0:110, diffssm, lty=2)
			lines(0:110,diffsc,col="#00DD11")
			lines(0:110,diffssmc,col="#00DD11",lty=2)
			lines(0:110,diffsa,col="#DD00DD")
			lines(0:110,diffssma,col="#DD00DD",lty=2)
			
			
		})
dev.off()
# ------------------------------------------------------------
# code used to decide which fitting and extrapolation cutoffs to use
#plotit <- function(X, a = 70,...){
#	ages <- a:110
#	X    <- X[(a + 1):nrow(X), ]
#	cols <- gray(seq(0,.7,length=5))
#	lwd  <- seq(1,3,length=5)
#	matplot(ages,X,type='l',col=cols,lty=1,lwd=lwd,log='y',...)
#}
#X <- Container.mx.m.75.80[,1,]
#yrs <- seq(1981,2011,by=10)
#for (i in 1:4){
#	plotit(Container.mx.f.70.85[,i,],main = paste0(yrs[i],"fit70, extrap 85"))
#	locator(1)
#	plotit(Container.mx.f.75.85[,i,],main =  paste0(yrs[i],"fit75, extrap 85"))
#	locator(1)
#}
#for (i in 1:4){
#	plotit(Container.mx.f.75.80[,i,],main = paste0(yrs[i],"fit75, extrap 80"))
#	locator(1)
#	plotit(Container.mx.f.75.85[,i,],main =  paste0(yrs[i],"fit75, extrap 85"))
#	locator(1)
#}
#plotit2 <- function(X, a = 70,...){
#	ages <- a:110
#	X    <- X[(a + 1):nrow(X), ]
#	cols <- gray(seq(0,.7,length=4))
#	lwd  <- seq(1,3,length=4)
#	matplot(ages,X,type='l',col=cols,lty=1,lwd=lwd,log='y',...)
#}
#for (i in 1:5){
#	plotit2(Container.mx.m.75.80[,,i],main = paste0(yrs[i],"fit75, extrap 80"))
#	locator(1)
#	plotit2(Container.mx.m.75.85[,,i],main =  paste0(yrs[i],"fit75, extrap 85"))
#	locator(1)
#}
#
#plot(75:84,SCO$N[SCO$age %in% 75:84 & SCO$sex == 1 & SCO$quintile == 1 & SCO$year == 1981], ylim = c(0,6000))
#lines(75:84,SCO$N[SCO$age %in% 75:84 & SCO$sex == 1 & SCO$quintile == 5 & SCO$year == 1981])
#
#ratio <- SCO$N[SCO$age %in% 75:84 & SCO$sex == 1 & SCO$quintile == 1 & SCO$year == 1991] / 
#		SCO$N[SCO$age %in% 75:84 & SCO$sex == 1 & SCO$quintile == 5 & SCO$year == 1991]
#
#plot(75:84,ratio)
#
#mx1 <- SCO$mx[SCO$sex == 1 & SCO$quintile == 1 & SCO$year == 1981] 
#mx2 <- SCO$mx[SCO$sex == 1 & SCO$quintile == 5 & SCO$year == 1981] 
#mx1sm <- SCO$mxsm[SCO$sex == 1 & SCO$quintile == 1 & SCO$year == 1981] 
#mx2sm <- SCO$mxsm[SCO$sex == 1 & SCO$quintile == 5 & SCO$year == 1981] 
#
#
#ratiomx <- mx1 / mx2
#ratiomxsm <- mx1sm / mx2sm
#plot(0:110,ratiomx)
#lines(0:110,ratiomxsm)
#for (i in 1:5){
#	plotit2(Container.mx.f.75.80[,,i],main = paste0(yrs[i],"fit75, extrap 80"))
#	locator(1)
#	plotit2(Container.mx.f.75.85[,,i],main =  paste0(yrs[i],"fit75, extrap 85"))
#	locator(1)
#}

# fitting using info from age 75-94 was the least screwy
# extrapolation after 85 rather than 80 seemed to preserve
# more natural variation.

#head(SCO)
#
#mx2e0cheap <- function(mx){
#	sum(exp(-cumsum(mx)))
#}
#
#mx2edcheap <- function(mx){
#	lx <- c(1,exp(-cumsum(mx)))
#	dx <- -diff(c(lx,0))
#	Lx <- (lx + c(lx[-1],0)) / 2
#	Tx <- rev(cumsum(rev(Lx)))
#	ex <- Tx / lx
#	sum(ex*dx)
#}
#SCO <- data.table(SCO)
#check <- SCO[,list(e0=mx2e0cheap(mx),e0sm=mx2e0cheap(mxsm)),
#		by=list(year,sex,quintile)]
#hist(check$e0-check$e0sm)
#check2 <- SCO[,list(ed0=mx2edcheap(mx),ed0sm=mx2edcheap(mxsm)),
#		by=list(year,sex,quintile)]
#
#hist(check2$ed0-check2$ed0sm)
#
#plot(SCO$mx[1:110],log='y')
head(SCO)