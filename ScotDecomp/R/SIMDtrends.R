# In this file we calculate all results for the full SIMD series
# this produces an estimate for each single year from 2001 to 2017 (17 years)
# Each step is repeated from both the full SIMD measure and the Income-only measure,
# for comparison. These results are not commensurable with the Carstairs series,
# and we judge that the main reason why the different deprivation measures produce
# different results is due to the modifiable areal unit problem: i.e. SIMD has 6000+
# geographic units, whereas Carstairs is tabulated on abotu 1/5 as many.


# header to determine working directory
me <- system("whoami",intern=TRUE)
if (me == "tim"){
	setwd("/home/tim/git/ScotDecomp/ScotDecomp")
}
if (me == "mpidr_d\\seaman"){
	setwd("U:/Conferences/PAA/2018 Denver/within and between/ScotDecomp")
	
}

# ---------------------

source("R/Functions.R")
library(foreign)
library(data.table)
# --------------------------

# read in data sets
# SIMD
SIMD <- read.dta("N:/Rosie/For Tim/BMJOpen/Data/Quintile_SIMD_Full.dta")
# Income-only
INC  <- read.dta("N:/Rosie/For Tim/BMJOpen/Data/Quintile_Income_Domain.dta")
head(SIMD)
# -----------------------------

# step 1: extrapolate to age 110
SIMD <- data.table(SIMD)
INC  <- data.table(INC)

# standardize names
setnames(SIMD,"SIMD_quintile","q")
setnames(INC,"income_quintile","q")

# extrapolate death rates to age 110 using modified HMD protocol
SIMD   <- SIMD[,extrap_dt(deaths_SIMD_quintile, population_SIMD_quintile),
		        by = list(yr, sex, q)]
	
INC    <- INC[,extrap_dt(deaths_income_quintile, population_income_quintile),
		        by = list(yr, sex, q)]
		
# standardize names	
setnames(SIMD,"X1","mx")
setnames(INC,"X1","mx")

# standardize names

# derive qx needed for matrix operations
SIMD$ax                       <- .5
INC$ax                        <- .5
# infants, same as SCO
SIMD$ax[SIMD$age == 0]        <- .1
INC$ax[INC$age == 0]          <- .1
# constant hazard close  out assumption
SIMD$ax[CARmx$age == 110]     <- 1 / SIMD$mx[CARmx$age == 110] 
INC$ax[INC$age == 110]        <- 1 / INC$mx[INC$age == 110] 

# now get qx:
SIMD$qx                       <- SIMD$mx / (1 + (1- SIMD$ax) * SIMD$mx )
INC$qx                        <- INC$mx / (1 + (1- INC$ax) * INC$mx )

# now do within between stuff

# generate B and W variance
SIMDBW <- SIMD[,list(age = unique(age),
				Bst = Byrsex(.SD), 
				Wst = Wyrsex(.SD)), 
		by = list(yr, sex)]

INCBW <- INC[,list(age = unique(age),
				Bst = Byrsex(.SD), 
				Wst = Wyrsex(.SD)), 
		by = list(yr, sex)]

# other quantities derives from B and B variance

# sum to total variance
SIMDBW$V       <- SIMDBW$B + SIMDBW$W
INCBW$V        <- INCBW$B + INCBW$W

# proportion between
SIMDBW$propB   <- SIMDBW$B / SIMDBW$V
INCBW$propB    <- INCBW$B / INCBW$V

# standard deviations
SIMDBW$sd      <- sqrt(SIMDBW$V)
INCBW$sd       <- sqrt(INCBW$V)

# combine results and save out
SIMDBW$measure <- "SIMD"
INCBW$measure  <- "INC"
BWout          <- rbind(SIMDBW, INCBW)
path           <- "N:/Rosie/For Tim/BMJOpen/Data/Results"
write.csv(BWout, file=file.path(path, "Between_SIMD_INC_Full.csv"))


# now sd against ex for quintiles, just age 0
e0sdSIMD  <- SIMD[, list(e0 = qxax2e0(qx, ax), sd = qx2sd0(qx)), by = list(yr, sex, q)]
e0sdINC   <- INC[, list(e0 = qxax2e0(qx, ax), sd = qx2sd0(qx)), by = list(yr, sex, q)]

path <-"N:/Rosie/For Tim/BMJOpen/Data/Results"
write.csv(e0sdSIMD, file=file.path(path,"EXSD_SIMD_Full.csv"))
write.csv(e0sdINC, file=file.path(path,"EXSD_INC_Full.csv"))


# exploratory plotting of results, deprecated

# quick look at age patterns over time:
#MPB <- acast(SIMDBW[SIMDBW$sex == "Male"], age~yr, value.var = "propB")
#FPB <- acast(SIMDBW[SIMDBW$sex == "Female"], age~yr, value.var = "propB")
#
#library(RColorBrewer)
#ramp <- colorRampPalette(brewer.pal(9,"PuRd"),space = "Lab")
#matplot(0:110,MPB,
#		type='l', xlim = c(0,85), lty = 1, col = ramp(17),ylim=c(0,.08),lwd=3)
#matplot(0:110,FPB,
#		type='l', xlim = c(0,85), lty = 1, col = ramp(17),ylim=c(0,.08),lwd=3)
#
#plot(2001:2017, MPB[1,],type='l')
#plot(2001:2017, MPB[10,],type='l')
#plot(2001:2017, MPB[30,],type='l')
#
#par(mfrow=c(1,2))
#matplot(2001:2017,t(MPB[as.character(seq(0,80,by=10)),]),
#		type='l', xlim = c(2001,2017), lty = 1, col = ramp(9),ylim=c(0,.08),lwd=3)
#matplot(2001:2017,t(FPB[as.character(seq(0,80,by=10)),]),
#		type='l', xlim = c(2001,2017), lty = 1, col = ramp(9),ylim=c(0,.08),lwd=3)

# again for income
#MPBI <- acast(INCBW[INCBW$sex == "Male"], age~yr, value.var = "propB")
#FPBI <- acast(INCBW[INCBW$sex == "Female"], age~yr, value.var = "propB")
#par(mfrow=c(1,2))
#matplot(0:110,MPBI,
#		type='l', xlim = c(0,85), lty = 1, col = ramp(17),ylim=c(0,.08),lwd=3)
#matplot(0:110,FPBI,
#		type='l', xlim = c(0,85), lty = 1, col = ramp(17),ylim=c(0,.08),lwd=3)
#
#

#display.brewer.all()
#cols <- brewer.pal(7,"YlOrBr")[3:7]
#colsf <- brewer.pal(7,"RdPu")[3:7]
#plot(NULL, type = "n", xlim = c(68,86), ylim = c(11,18), ann = FALSE, asp=2,las=1)
#for(i in 1:5){
#	ind <- e0sdSIMD$sex == "Male" & e0sdSIMD$q == i
#	lines(e0sdSIMD$e0[ind],e0sdSIMD$sd[ind],lwd=2,col = cols[i])
#	
#	ind <- e0sdSIMD$sex == "Female" & e0sdSIMD$q == i
#	lines(e0sdSIMD$e0[ind],e0sdSIMD$sd[ind],lwd=2,col = colsf[i])
#	
#}
#
#
#
#
#desc <- read.csv("Data/descriptives.csv",stringsAsFactors = FALSE)
#head(desc)
#
## little spot checks
#
#
#
#plot(NULL, type= "n", xlim = c(65,85),ylim=c(10,20))
#for (q in 1:5){
#ind1 <- desc$sex == 2 & desc$quintile_2 == q & desc$age == 0 & desc$year %in% c(2001, 2011)
#lines(desc$ex[ind1],desc$sd[ind1],col = "red", lwd=2)
#points(desc$ex[ind1],desc$sd[ind1],col = "red",pch=16)
#
#ind2 <- e0sdSIMD$sex == "Female" & e0sdSIMD$q == 6-q & e0sdSIMD$yr %in% c(2001, 2011)
#lines(e0sdSIMD$e0[ind2],e0sdSIMD$sd[ind2],col = "blue", lwd=2)
#points(e0sdSIMD$e0[ind2],e0sdSIMD$sd[ind2],col = "blue",pch=16)
#
#ind2 <- e0sdINC$sex == "Female" & e0sdINC$q == 6-q & e0sdINC$yr %in% c(2001, 2011)
#lines(e0sdINC$e0[ind2],e0sdINC$sd[ind2],col = "green", lwd=2)
#points(e0sdINC$e0[ind2],e0sdINC$sd[ind2],col = "green",pch=16)
#}
#
#
#
## time trends in e0:
#
#plot(NULL, type= "n", xlim = c(1981, 2017), ylim = c(65,87),las=1,ylab="e(0)",xlab="Year")
#for (q in c(1,3,5)){
#ind1 <- desc$sex == 1 & desc$quintile_2 == q & desc$age == 0
#lines(c(1981,1991,2001,2011),desc$ex[ind1],col = "red",lwd=2 )
#ind2 <- e0sdSIMD$sex == "Male" & e0sdSIMD$q == 6-q 
#lines(2001:2017,e0sdSIMD$e0[ind2],col = "red",lty="8282",lwd=2)
#ind3 <- e0sdINC$sex == "Male" & e0sdINC$q == 6-q 
#lines(2001:2017,e0sdINC$e0[ind3],col = "red",lty="4242",lwd=2)
#segments(2011,e0c$car[e0c$sex == "Male" & e0c$q %in% c(1,3,5)],
#		2011, e0c$simd1[e0c$sex == "Male" & e0c$q %in% c(1,3,5)])
#}
## time trends in sd
#
#plot(NULL, type= "n", xlim = c(1981, 2017), ylim = c(10,20),las=1,ylab="e(0)",xlab="Year")
#for (q in c(1,3,5)){
#ind1 <- desc$sex == 1 & desc$quintile_2 == q & desc$age == 0
#lines(c(1981,1991,2001,2011),desc$sd[ind1],col = "red",lwd=2 )
#ind2 <- e0sdSIMD$sex == "Male" & e0sdSIMD$q == 6-q 
#lines(2001:2017,e0sdSIMD$sd[ind2],col = "red",lty="8282",lwd=2)
##ind3 <- e0sdINC$sex == "Male" & e0sdINC$q == 6-q 
##lines(2001:2017,e0sdINC$sd[ind3],col = "red",lty="4242",lwd=2)
#
#}
#segments(2011,e0c$car[e0c$sex == "Male" & e0c$q %in% c(1,3,5)],
#		2011, e0c$simd1[e0c$sex == "Male" & e0c$q %in% c(1,3,5)])
#
