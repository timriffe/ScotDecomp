
# Author: tim
###############################################################################
# a script to extrapolate mx from age 85+ to age 110+

source(file.path("R","0_Functions.R"))

# read in data
SCO <- read.csv(file.path("Data","Inputs","lifetables_deciles_scotland.csv"))

# account for different year groupings
SCO$Fac                   <- 3
SCO$Fac[SCO$year == 1991] <- 2

# in case of sims, we want to be on the right scale.
SCO$N <- SCO$N / SCO$Fac
SCO$D <- SCO$D / SCO$Fac
# get maxlik kannisto functions
SCO <- SCO[SCO$decile != 999, ]
ind <- SCO$year == 2011
SCO$decile[ind] <- abs(SCO$decile[ind] - 11)
#matplot(acast(SCO[SCO$sex == 1 & SCO$age == 0, ], decile~year, value.var ="ex"),type='l')
# this makes arrays, each layer is a quintile consisting in an AP matrix
DXm   <- acast(SCO[SCO$sex == 1, ], age~year~decile, value.var = "D")
DXf   <- acast(SCO[SCO$sex == 2, ], age~year~decile, value.var = "D")
NXm   <- acast(SCO[SCO$sex == 1, ], age~year~decile, value.var = "N")
NXf   <- acast(SCO[SCO$sex == 2, ], age~year~decile, value.var = "N")

# boxes to put results in
Container.mx.m <- array(dim=c(111,dim(DXm)[2:3]),dimnames=list(0:110,seq(1981,2011,10),c(1:10)))
Container.mx.f <- array(dim=c(111,dim(DXf)[2:3]),dimnames=list(0:110,seq(1981,2011,10),c(1:10)))

# extrapolation loop.
# just ignore warnings and printouts
for (i in 1:dim(DXm)[3]){
	Container.mx.m[,,i] <- extrap(DXm[,,i],NXm[,,i])
	Container.mx.f[,,i] <- extrap(DXf[,,i],NXf[,,i])
}

# reshape to long for easier operations
mxmlong     <- melt(Container.mx.m, value.name="mx",varnames=c("age","year","decile"))
mxflong     <- melt(Container.mx.f, value.name="mx",varnames=c("age","year","decile"))

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
#save(SCOlong, file = file.path("Data","Derived","SCOlongDeciles.Rdata"))
Byrsex <- function(.SD,w="stationary"){
	QXkts <- acast(.SD, age ~decile, value.var = "qx")
	if (w == "stationary"){
		out <- Vbetween(QXkts)
	}
	if (w == "raw"){
		NXkts <- acast(.SD, age ~decile, value.var = "N")
		raww  <- NXkts / rowSums(NXkts)
		out   <- Vbetween(QXkts,raww)
	}
	out
}
Wyrsex <- function(.SD,w="stationary"){
	QXkts <- acast(.SD, age ~ decile, value.var = "qx")
	if (w == "stationary"){
		out <- c(Vwithin(QXkts))
	}
	if (w == "raw"){
		NXkts <- acast(.SD, age ~ decile, value.var = "N")
		raww  <- NXkts / rowSums(NXkts)
		out <- c(Vwithin(QXkts,raww))
	}
	out
}
# convert data type for easier calcs
SCO <- data.table(SCOlong)

# calculate between and within components directly. By year and sex,
# implying 4 x 2 subsets.
SCOB <- SCO[,list(age = unique(age),
				Bst = Byrsex(.SD), 
				Wst = Wyrsex(.SD)), 
		by = list(year, sex)]

# total is the sum
SCOB$V     <- SCOB$B + SCOB$W

# proportion between (reported in manuscript)
SCOB$propB <- SCOB$B / SCOB$V

# standard deviation (reported in manuscript)
SCOB$sd    <- sqrt(SCOB$V)

# convert back to standard data.frame, just because
SCOB       <- as.data.frame(SCOB)

# save out intermediate data object for downstream analyses, viz
write.csv(SCOB, file = file.path("Data","Derived","DecileCompare.csv"))

# quick look
do.this <- FALSE
if (do.this){
X <- acast(SCOB, age ~ year ~ sex, value.var = "propB")

matplot(X[,,2], type = 'l',ylim=c(0,.04))
}

# end