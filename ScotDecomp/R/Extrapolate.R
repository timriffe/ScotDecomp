
# Author: tim
###############################################################################

# a script to extrapolate mx from age 85+ to age 110+
# header to determine working directory
me <- system("whoami")
if (me == "tim"){
	setwd("/home/tim/git/ScotDecomp/ScotDecomp")
}
if (me == "mpidr_d\seaman"){
	setwd("U:\Conferences\PAA\2018 Denver\within and between\ScotDecomp")
}

# ---------------------

source("R/Functions.R")

# read in data
SCO <- read.csv("Data/lifetables_quintiles_scotland.csv")
SCO$Fac <- 3
SCO$Fac[SCO$year == 1991] <- 2

SCO$N <- SCO$N / SCO$Fac
SCO$D <- SCO$D / SCO$Fac
# get maxlik kannisto functions


library(reshape2)

# this makes arrays, each layer is a quintile consisting in an AP matrix
DXm <- acast(SCO[SCO$sex == 1, ], age~year~quintile_2, value.var = "D")
DXf <- acast(SCO[SCO$sex == 2, ], age~year~quintile_2, value.var = "D")
NXm <- acast(SCO[SCO$sex == 1, ], age~year~quintile_2, value.var = "N")
NXf <- acast(SCO[SCO$sex == 2, ], age~year~quintile_2, value.var = "N")


# a function to deal with a layer of an array
extrap <- function(DX, NX){
	DX <- DX[-nrow(DX), ]
	NX <- NX[-nrow(NX), ]
	mat110 <- matrix(NA,ncol=ncol(DX),nrow=111,dimnames=list(0:110,colnames(DX)))
	DXe <- mat110
	DXe[rownames(DX),colnames(DX)] <- DX
	NXe <- mat110
	NXe[rownames(NX),colnames(NX)] <- NX
	
	ages.i <- rep(nrow(DX),ncol(DX)) + 1
	ltper_mx_v5(DXe,NXe,ages.i,fit.from.i = 76 )
}

Container.mx.m <- array(dim=c(111,dim(DXm)[2:3]),dimnames=list(0:110,seq(1981,2011,10),c(1:5,999)))
Container.mx.f <- array(dim=c(111,dim(DXf)[2:3]),dimnames=list(0:110,seq(1981,2011,10),c(1:5,999)))

# just ignore warnings and printouts
for (i in 1:dim(DXm)[3]){
	Container.mx.m[,,i] <- extrap(DXm[,,i],NXm[,,i])
	Container.mx.f[,,i] <- extrap(DXf[,,i],NXf[,,i])
}
# reshape to long for easier operations
mxmlong     <- melt(Container.mx.m, value.name="mx",varnames=c("age","year","quintile_2"))
mxflong     <- melt(Container.mx.f, value.name="mx",varnames=c("age","year","quintile_2"))
mxmlong$sex <- 1
mxflong$sex <- 2
SCOlong     <- rbind(mxmlong, mxflong)

# now get back to qx
SCOlong$ax                     <- .5
# infants, same as SCO
SCOlong$ax[SCOlong$age == 0]   <- .1
# constant hazard closeout assumption
SCOlong$ax[SCOlong$age == 110] <- 1 / SCOlong$mx[SCOlong$age == 110] 

# now get qx:
SCOlong$qx <- SCOlong$mx / (1 + (1- SCOlong$ax) * SCOlong$mx )
# closeout:
SCOlong$qx[SCOlong$age == 110] <- 1

# quick check to see if anything stands out:
# matplot(0:110,acast(SCOlong, age~year+quintile_2+sex,value.var="mx"),log='y',type='l')

save(SCOlong, file = "Data/SCOlong.Rdata")
# end