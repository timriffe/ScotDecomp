
# this is 2011, note
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
data.path <- "N:/Rosie/For Tim/BMJOpen/Data/V12"
Car   <- read.dta(file.path(data.path,"Carstairs_score_for_datazones_V12.dta"))
Car   <- data.table(Car)
SIMD1 <- Car[,list(D = sum(total_death), E = sum(population)), by = list(sex,ageyrs,simd2012_sc_quintile)]
SIMD2 <- Car[,list(D = sum(total_death), E = sum(population)), by = list(sex,ageyrs,quintile_income)]
CAR   <- Car[,list(D = sum(total_death), E = sum(population)), by = list(sex,ageyrs,quintile_carstairs_for_datazone)]

setnames(SIMD1,"simd2012_sc_quintile","q")
setnames(SIMD2,"quintile_income","q")
setnames(CAR,"quintile_carstairs_for_datazone","q")
setnames(SIMD1,"ageyrs","age")
setnames(SIMD2,"ageyrs","age")
setnames(CAR,"ageyrs","age")

CAR   <- CAR[!is.na(q)]
CAR   <- CAR[order(sex,q,age)]
SIMD1 <- SIMD1[order(sex,q,age)]
SIMD2 <- SIMD2[order(sex,q,age)]


extrap_dt <- function(D,E){
	
	D      <- as.matrix(D[1:90])
	E      <- as.matrix(E[1:90])
	rownames(D) <- 0:89
	rownames(E) <- 0:89
	colnames(D) <- 1
	colnames(E) <- 1
	agenew <- 90:110
	extend <- matrix(NA,nrow=length(agenew),ncol=ncol(D),dimnames=list(agenew,colnames(D)))
	D      <- rbind(D,extend)
	E      <- rbind(E,extend)
	Mx     <- ltper_mx_v5(D,E, fit.from.i = 76)
	data.frame(age = 0:110,mx=Mx)
}

CARmx   <- CAR[,extrap_dt(D,E),by=list(sex,q)]
SIMD1mx <- SIMD1[,extrap_dt(D,E),by=list(sex,q)]
SIMD2mx <- SIMD2[,extrap_dt(D,E),by=list(sex,q)]

setnames(CARmx,"X1","mx")
setnames(SIMD1mx,"X1","mx")
setnames(SIMD2mx,"X1","mx")
#plot(0:110, CARmx$X1[CARmx$q==3 & CARmx$sex=="Female"],log="y",col="black",type='l')
#lines(0:110, SIMD1mx$X1[SIMD1mx$q==3 & SIMD1mx$sex=="Female"],col="blue")
#lines(0:110, SIMD2mx$X1[SIMD2mx$q==3 & SIMD2mx$sex=="Female"],col="red")
#
#lines(0:110, CARmx$X1[CARmx$q==1 & CARmx$sex=="Female"],col="black",lty=2)
#lines(0:110, SIMD1mx$X1[SIMD1mx$q==1 & SIMD1mx$sex=="Female"],col="blue",lty=2)
#lines(0:110, SIMD2mx$X1[SIMD2mx$q==1 & SIMD2mx$sex=="Female"],col="red",lty=2)
#
#lines(0:110, CARmx$X1[CARmx$q==5 & CARmx$sex=="Female"],col="black",lty="8282")
#lines(0:110, SIMD1mx$X1[SIMD1mx$q==5 & SIMD1mx$sex=="Female"],col="blue",lty="8282")
#lines(0:110, SIMD2mx$X1[SIMD2mx$q==5 & SIMD2mx$sex=="Female"],col="red",lty="8282")


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
# Modified from 2_Analysis.R
# redo analysis for 2011, all 3 quantile groupings
Byrsex <- function(.SD,w="stationary"){
	QXkts <- acast(.SD, age ~ q, value.var = "qx")
	out <- Vbetween(QXkts)
	out
}
Wyrsex <- function(.SD){
	QXkts <- acast(.SD, age ~q, value.var = "qx")
	out  <- c(Vwithin(QXkts))
	out
}

# generate B and W variance
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
# sum to total variance
CARBW$V       <- CARBW$B + CARBW$W
SIMD1BW$V     <- SIMD1BW$B + SIMD1BW$W
SIMD2BW$V     <- SIMD2BW$B + SIMD2BW$W
# proportion between
CARBW$propB   <- CARBW$B / CARBW$V
SIMD1BW$propB <- SIMD1BW$B / SIMD1BW$V
SIMD2BW$propB <- SIMD2BW$B / SIMD2BW$V
# standard deviations
CARBW$sd      <- sqrt(CARBW$V)
SIMD1BW$sd    <- sqrt(SIMD1BW$V)
SIMD2BW$sd    <- sqrt(SIMD2BW$V)

# something easier to plot?
CARp         <- acast(CARBW, age~sex, value.var = "propB")
SIMD1p         <- acast(SIMD1BW, age~sex, value.var = "propB")
SIMD2p         <- acast(SIMD2BW, age~sex, value.var = "propB")

a <- 0:110

# SI figure
pdf("Figures/BetweenPropCARSIMD.pdf")
plot(a, CARp[,2], type = 'l', col = "black", xlim = c(0,85), 
		ylim = c(0,.07),las=1,main="Proportion of variance due to differences between quintiles in 2011
Carstairs and SIMD",
xlab = "Age",ylab = "Proportion")
lines(a, SIMD1p[,2], col = "blue")
lines(a, SIMD2p[,2], col = "red")

lines(a, CARp[,1], col = "black",lty=5)
lines(a, SIMD1p[,1], col = "blue",lty=5)
lines(a, SIMD2p[,1], col = "red",lty=5)

legend("topright",
		lty=c(1,1,1,5,5,5),
		col=c("black","blue","red","black","blue","red"),
		legend = c("Males Carstairs","Males SIMD","Males SIMD(inc)",
				"Females Carstairs","Females SIMD","Females SIMD(inc)"),
		bty="n"
)
dev.off()

# ----------------------------------------------------
# age 0 SD for each
rbind(unlist(CARBW[CARBW$age==0,"sd"]),
		unlist(SIMD1BW[SIMD1BW$age==0,"sd"]),
		unlist(SIMD2BW[SIMD2BW$age==0,"sd"]))

# for curiosity, what is e0 of each quintile
e0c  <- CARmx[,list(car=qxax2e0(qx,ax)),by=list(sex,q)]
e0s1 <- SIMD1mx[,list(simd1=qxax2e0(qx,ax)),by=list(sex,q)]
e0s2 <- SIMD2mx[,list(simd2=qxax2e0(qx,ax)),by=list(sex,q)]

e0c$simd1 <- e0s1$simd1
e0c$simd2 <- e0s2$simd2

# figure emailed to Rosie 19-10-2018
png("Figures/SIMDe0compare.png")
plot(e0c$q,e0c$car,pch=1,col=rep(c("red","blue"),each=5))
points(e0c$q,e0c$simd1,pch=2,col=rep(c("red","blue"),each=5))
points(e0c$q,e0c$simd2,pch=3,col=rep(c("red","blue"),each=5))
legend("bottomright",col=rep(c("red","blue"),each=3),pch=c(1,2,3,1,2,3),
		legend=c("Males Carstairs","Males SIMD","Males SIMD(inc)",
				"Females Carstairs","Females SIMD","Females SIMD(inc)"))
dev.off()


