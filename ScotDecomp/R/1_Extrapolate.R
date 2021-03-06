
# Author: tim
###############################################################################
# a script to extrapolate mx from age 85+ to age 110+

source(file.path("R","0_Functions.R"))

# read in data
SCO <- read.csv(file.path("Data","Inputs","lifetables_quintiles_scotland.csv"))

# account for different year groupings
SCO$Fac                   <- 3
SCO$Fac[SCO$year == 1991] <- 2

# in case of sims, we want to be on the right scale.
SCO$N <- SCO$N / SCO$Fac
SCO$D <- SCO$D / SCO$Fac
# get maxlik kannisto functions

# this makes arrays, each layer is a quintile consisting in an AP matrix
DXm   <- acast(SCO[SCO$sex == 1, ], age~year~quintile_2, value.var = "D")
DXf   <- acast(SCO[SCO$sex == 2, ], age~year~quintile_2, value.var = "D")
NXm   <- acast(SCO[SCO$sex == 1, ], age~year~quintile_2, value.var = "N")
NXf   <- acast(SCO[SCO$sex == 2, ], age~year~quintile_2, value.var = "N")

# boxes to put results in
Container.mx.m <- array(dim=c(111,dim(DXm)[2:3]),dimnames=list(0:110,seq(1981,2011,10),c(1:5,999)))
Container.mx.f <- array(dim=c(111,dim(DXf)[2:3]),dimnames=list(0:110,seq(1981,2011,10),c(1:5,999)))

# extrapolation loop.
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
#matplot(0:110,acast(SCOlong, age~year+quintile_2+sex,value.var="mx"),log='y',type='l')

# save out for downstream calcs
save(SCOlong, file = file.path("Data","Derived","SCOlong.Rdata"))

# end