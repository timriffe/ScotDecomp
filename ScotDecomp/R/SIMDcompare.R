
# Here we have death and population counts in SIMD datazones, 
# but with separate population-weighted quintile classification
# by Carstairs, SIMD, and the Income-only domain of SIMD. This
# allows for a measure comparison without free of the modifiable
# areal unit problem.


# 1) deaths grouped for years 2010,2011,2012
# 2) census in 2011, counts multiplied by 3

me <- system("whoami",intern=TRUE)
if (me == "tim"){
	setwd("/home/tim/git/ScotDecomp/ScotDecomp")
} 
if (me == "mpidr_d\\riffe"){
	setwd("U:/git/ScotDecomp/ScotDecomp")
}
source("R/Functions.R")
library(foreign)
library(data.table)
library(reshape2)

# Data kept on N drive, MPIDR PC only.
data.path <- "N:/Rosie/For Tim/BMJOpen/Data/V12"
Car       <- read.dta(file.path(data.path,"Carstairs_score_for_datazones_V12.dta"))
Car       <- data.table(Car)

# aggregate within the quintiles of each measure; creates 3 datasets.
# from now on each coding step repeated 3 times, sorry!
SIMD1 <- Car[,list(D = sum(total_death), E = sum(population)), by = list(sex, ageyrs, simd2012_sc_quintile)]
SIMD2 <- Car[,list(D = sum(total_death), E = sum(population)), by = list(sex, ageyrs, quintile_income)]
CAR   <- Car[,list(D = sum(total_death), E = sum(population)), by = list(sex, ageyrs, quintile_carstairs_for_datazone)]

# standardize names
setnames(SIMD1,"simd2012_sc_quintile","q")
setnames(SIMD2,"quintile_income","q")
setnames(CAR,"quintile_carstairs_for_datazone","q")
setnames(SIMD1,"ageyrs","age")
setnames(SIMD2,"ageyrs","age")
setnames(CAR,"ageyrs","age")

# some Carstairs quintiles were NA, remove
CAR     <- CAR[!is.na(q)]

# order each dataset
CAR     <- CAR[order(sex, q, age)]
SIMD1   <- SIMD1[order(sex, q, age)]
SIMD2   <- SIMD2[order(sex, q, age)]

# extrapolate mx to age 110 (each sex,quintile)
CARmx   <- CAR[, extrap_dt(D, E), by = list(sex,q)]
SIMD1mx <- SIMD1[, extrap_dt(D, E), by = list(sex,q)]
SIMD2mx <- SIMD2[, extrap_dt(D, E), by = list(sex,q)]

# standardize names
setnames(CARmx,"X1","mx")
setnames(SIMD1mx,"X1","mx")
setnames(SIMD2mx,"X1","mx")

# now get back to qx
CARmx$ax                       <- .5
SIMD1mx$ax                     <- .5
SIMD2mx$ax                     <- .5

# infants, same as SCO
CARmx$ax[CARmx$age == 0]       <- .1
SIMD1mx$ax[SIMD1mx$age == 0]   <- .1
SIMD2mx$ax[SIMD2mx$age == 0]   <- .1

# constant hazard closeout assumption
CARmx$ax[CARmx$age == 110]     <- 1 / CARmx$mx[CARmx$age == 110] 
SIMD1mx$ax[SIMD1mx$age == 110] <- 1 / SIMD1mx$mx[SIMD1mx$age == 110] 
SIMD2mx$ax[SIMD2mx$age == 110] <- 1 / SIMD2mx$mx[SIMD2mx$age == 110] 

# now get qx:
CARmx$qx                       <- CARmx$mx / (1 + (1- CARmx$ax) * CARmx$mx )
SIMD1mx$qx                     <- SIMD1mx$mx / (1 + (1- SIMD1mx$ax) * SIMD1mx$mx )
SIMD2mx$qx                     <- SIMD2mx$mx / (1 + (1- SIMD2mx$ax) * SIMD2mx$mx )

# closeout:
CARmx$qx[CARmx$age == 110]     <- 1
SIMD1mx$qx[SIMD1mx$age == 110] <- 1
SIMD2mx$qx[SIMD2mx$age == 110] <- 1

# now we have all the pieces we need for analysis

# generate Between and Within variance
CARBW <- CARmx[,list(age = unique(age),
				Bst = Byrsex(.SD), 
				Wst = Wyrsex(.SD)), 
		       by = list(sex)]
SIMD1BW <- SIMD1mx[,list(age = unique(age),
					   Bst = Byrsex(.SD), 
					   Wst = Wyrsex(.SD)), 
			   by = list(sex)]
SIMD2BW <- SIMD2mx[,list(age = unique(age),
					   Bst = Byrsex(.SD), 
					   Wst = Wyrsex(.SD)), 
			   by = list(sex)]
	   
# get measures detived from between and within variance
# sum to total variance
CARBW$V         <- CARBW$B + CARBW$W
SIMD1BW$V       <- SIMD1BW$B + SIMD1BW$W
SIMD2BW$V       <- SIMD2BW$B + SIMD2BW$W
# proportion between
CARBW$propB     <- CARBW$B / CARBW$V
SIMD1BW$propB   <- SIMD1BW$B / SIMD1BW$V
SIMD2BW$propB   <- SIMD2BW$B / SIMD2BW$V
# standard deviations
CARBW$sd        <- sqrt(CARBW$V)
SIMD1BW$sd      <- sqrt(SIMD1BW$V)
SIMD2BW$sd      <- sqrt(SIMD2BW$V)

# add new column to stack results
CARBW$measure   <- "Carstairs"
SIMD1BW$measure <- "SIMD"
SIMD2BW$measure <- "INC"

CompareOut <- rbind(CARBW, SIMD1BW, SIMD2BW)
path <-"N:/Rosie/For Tim/BMJOpen/Data/Results"
write.csv(CompareOut, file = file.path(path, "Between_Compare_2011.csv"))

# Derive e0 and sd0 for each quintile also to report
e0c          <- CARmx[, list(e0 = qxax2e0(qx, ax), sd = qx2sd0(qx)), by = list(sex, q)]
e0s1         <- SIMD1mx[, list(e0 = qxax2e0(qx, ax), sd = qx2sd0(qx)), by = list(sex, q)]
e0s2         <- SIMD2mx[, list(e0 = qxax2e0(qx, ax), sd = qx2sd0(qx)), by = list(sex, q)]

e0c$measure  <- "Carstairs"
e0s1$measure <- "SIMD"
e0s2$measure <- "INC"

e0out        <- rbind(e0c,e0s1, e0s2)

path         <-"N:/Rosie/For Tim/BMJOpen/Data/Results"
write.csv(e0out,file=file.path(path,"EXSD_Compare_2011.csv"))

# exploratory plots and other assorted deprecated code
# ----------------------------------------------------
#
# something easier to plot?
#CARp         <- acast(CARBW, age~sex, value.var = "propB")
#SIMD1p       <- acast(SIMD1BW, age~sex, value.var = "propB")
#SIMD2p       <- acast(SIMD2BW, age~sex, value.var = "propB")

#a <- 0:110

# SI figure
#pdf("Figures/BetweenPropCARSIMD.pdf")
#plot(a, CARp[,2], type = 'l', col = "black", xlim = c(0,85), 
#		ylim = c(0,.07),las=1,main="Proportion of variance due to differences between quintiles in 2011
#Carstairs and SIMD",
#xlab = "Age",ylab = "Proportion")
#lines(a, SIMD1p[,2], col = "blue")
#lines(a, SIMD2p[,2], col = "red")
#
#lines(a, CARp[,1], col = "black",lty=5)
#lines(a, SIMD1p[,1], col = "blue",lty=5)
#lines(a, SIMD2p[,1], col = "red",lty=5)
#
#legend("topright",
#		lty=c(1,1,1,5,5,5),
#		col=c("black","blue","red","black","blue","red"),
#		legend = c("Males Carstairs","Males SIMD","Males SIMD(inc)",
#				"Females Carstairs","Females SIMD","Females SIMD(inc)"),
#		bty="n"
#)
#dev.off()

# ----------------------------------------------------
# age 0 SD for each
#rbind(unlist(CARBW[CARBW$age==0,"sd"]),
#		unlist(SIMD1BW[SIMD1BW$age==0,"sd"]),
#		unlist(SIMD2BW[SIMD2BW$age==0,"sd"]))


#
#e0c$simd1     <- e0s1$simd1
#e0c$simd2     <- e0s2$simd2
#
#
#
## figure emailed to Rosie 19-10-2018
#png("Figures/SIMDe0compare.png")
#plot(e0c$q,e0c$car,pch=1,col=rep(c("red","blue"),each=5))
#points(e0c$q,e0c$simd1,pch=2,col=rep(c("red","blue"),each=5))
#points(e0c$q,e0c$simd2,pch=3,col=rep(c("red","blue"),each=5))
#legend("bottomright",col=rep(c("red","blue"),each=3),pch=c(1,2,3,1,2,3),
#		legend=c("Males Carstairs","Males SIMD","Males SIMD(inc)",
#				"Females Carstairs","Females SIMD","Females SIMD(inc)"))
#dev.off()

# -------------------------------------
# This was a once-off spot check to make sure that the counts
# were commemnsurable between the 3 derived datasets. They are:

do.this <- FALSE
if (do.this){
	SIMD1T <- Car[,list(D = sum(total_death), E = sum(population)), by = list(sex,ageyrs)]
	SIMD2T <- Car[,list(D = sum(total_death), E = sum(population)), by = list(sex,ageyrs)]
	CART   <- Car[,list(D = sum(total_death), E = sum(population)), by = list(sex,ageyrs)]
	
	CART   <- CART[order(sex, ageyrs)]
	SIMD1T <- SIMD1T[order(sex, ageyrs)]
	SIMD2T <- SIMD2T[order(sex, ageyrs)]
	
    # extrap
	CART   <- CART[, extrap_dt(D, E), by = list(sex)]
	SIMD1T <- SIMD1T[, extrap_dt(D, E), by = list(sex)]
	SIMD2T <- SIMD2T[, extrap_dt(D, E), by= list(sex)]
	
	setnames(CART,"X1","mx")
	setnames(SIMD1T,"X1","mx")
	setnames(SIMD2T,"X1","mx")
	
    # now get back to qx
	CART$ax                       <- .5
	SIMD1T$ax                     <- .5
	SIMD2T$ax                     <- .5
    # infants, same as SCO
	CART$ax[CART$age == 0]       <- .1
	SIMD1T$ax[SIMD1T$age == 0]   <- .1
	SIMD2T$ax[SIMD2T$age == 0]   <- .1
    # constant hazard closeout assumption
	CART$ax[CART$age == 110]     <- 1 / CART$mx[CART$age == 110] 
	SIMD1T$ax[SIMD1T$age == 110] <- 1 / SIMD1T$mx[SIMD1T$age == 110] 
	SIMD2T$ax[SIMD2T$age == 110] <- 1 / SIMD2T$mx[SIMD2T$age == 110] 
	
    # now get qx:
	CART$qx                       <- CART$mx / (1 + (1- CART$ax) * CART$mx )
	SIMD1T$qx                     <- SIMD1T$mx / (1 + (1- SIMD1T$ax) * SIMD1T$mx )
	SIMD2T$qx                     <- SIMD2T$mx / (1 + (1- SIMD2T$ax) * SIMD2T$mx )
    # closeout:
	CART$qx[CART$age == 110]     <- 1
	SIMD1T$qx[SIMD1T$age == 110] <- 1
	
	# now life expectancy results
	CART[, list(car = qxax2e0(qx, ax)), by = list(sex)]
	SIMD1T[, list(simd1 = qxax2e0(qx, ax)), by = list(sex)]
	SIMD2T[, list(simd2 = qxax2e0(qx, ax)), by = list(sex)]
}

# end
