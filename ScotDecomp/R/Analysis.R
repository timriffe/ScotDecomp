# Author: tim
###############################################################################

# TODO: make graphs of within and between variance


# header to determine working directory
me <- system("whoami",intern=TRUE)
if (me == "tim"){
	setwd("/home/tim/git/ScotDecomp/ScotDecomp")
}
if (me == "mpidr_d\\seaman"){
	setwd("U:/Conferences/PAA/2018 Denver/within and between/ScotDecomp")
	
}

# ---------------------
library(reshape2)
library(data.table)
source("R/Functions.R")
#source("R/Extrapolate.R")


SCO <- local(get(load("Data/SCOlong.Rdata")))
SCO <- SCO[SCO$quintile_2 != 999, ]

# calculate things manually (too much typing sorry)

# make a function that calcs within and between for a year and sex:
# raw weighting not prepared because N in orig dat aonly up to 85+,
# could resistribute 85+ according to l(85+), but that's for another
# day. 
Byrsex <- function(.SD,w="stationary"){
	QXkts <- acast(.SD, age ~quintile_2, value.var = "qx")
	if (w == "stationary"){
		out <- Vbetween(QXkts)
	}
	if (w == "raw"){
		NXkts <- acast(.SD, age ~quintile_2, value.var = "N")
		raww  <- NXkts / rowSums(NXkts)
		out <- Vbetween(QXkts,raww)
	}
	out
}
Wyrsex <- function(.SD,w="stationary"){
	QXkts <- acast(.SD, age ~quintile_2, value.var = "qx")
	if (w == "stationary"){
		out <- c(Vwithin(QXkts))
	}
	if (w == "raw"){
		NXkts <- acast(.SD, age ~quintile_2, value.var = "N")
		raww  <- NXkts / rowSums(NXkts)
		out <- c(Vwithin(QXkts,raww))
	}
	out
}

SCO <- data.table(SCO)

SCOB <- SCO[,list(age = unique(age),
				  Bst = Byrsex(.SD), 
				  Wst = Wyrsex(.SD)), 
		  by = list(year, sex)]

SCOB$V     <- SCOB$B + SCOB$W
SCOB$propB <- SCOB$B / SCOB$V
SCOB$sd    <- sqrt(SCOB$V)

SCOB       <- as.data.frame(SCOB)

save(SCOB,file="Data/SCOB.Rdata")


mp         <- acast(SCOB[SCOB$sex == 1, ], age~year, value.var = "propB")
fp         <- acast(SCOB[SCOB$sex == 2, ], age~year, value.var = "propB")
mp         <- mp[1:86, ]
fp         <- fp[1:86, ]
# from here down not reworked yet
a          <- 0:85

graphics.off()
pdf("Figures/BetweenPropMales.pdf")
matplot(a, mp, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "Proportion",xlab = "Age", ylim=c(0,.04),
		main = "",
		las = 1,
		cex.lab = 1.4)
abline(v=35)
text(20,mp[21, ], c(1981,1991,2001,2011),pos=3,cex=1)
dev.off()

pdf("Figures/BetweenPropFemales.pdf")
matplot(a, fp, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "Proportion",xlab = "Age", ylim=c(0,.04),
		main = "",
		las = 1,
		cex.lab = 1.4)
abline(v=35)
text(40,fp[41, ], c(1981,1991,2001,2011),pos=c(1,3,3,3),cex=1)
dev.off()


# ----------------------------
# sd:

msd         <- acast(SCOB[SCOB$sex == 1, ], age~year, value.var = "sd")
fsd         <- acast(SCOB[SCOB$sex == 2, ], age~year, value.var = "sd")
msd         <- msd[1:86, ]
fsd         <- fsd[1:86, ]

pdf("Figures/TotalsdMales.pdf")
matplot(a, msd, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "sd",xlab = "Age", ylim = c(0,16),
		main = "",
		las = 1,
		cex.lab = 1.4)
text(20,msd[21, "1981" ],1981,pos=2) 
text(30,msd[31, "2011" ],2011,pos=4) 
dev.off()

pdf("Figures/TotalsdFemales.pdf")
matplot(a, fsd, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "sd",xlab = "Age", ylim=c(0,16),
		main = "",
		las = 1,
		cex.lab = 1.4)

text(70,fsd[71, "1981"], 1981, pos=2)
text(70,fsd[71, "2011"], 2011, pos=4)
dev.off()

# new graphs as of 17-Nov-2017
# absolute var between:
mb         <- acast(SCOB[SCOB$sex == 1, ], age~year, value.var = "Bst")
fb         <- acast(SCOB[SCOB$sex == 2, ], age~year, value.var = "Bst")
mb         <- mb[1:86, ]
fb         <- fb[1:86, ]

pdf("Figures/TotalbstMales.pdf")
matplot(a, mb, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "between variance",xlab = "Age", ylim = c(0,9),
		main = "",
		las = 1,
		cex.lab = 1.4)
text(20,mb[21,]-.3,c(1981,1991,2001,2011)) 
dev.off()

pdf("Figures/TotalbstFemales.pdf")
matplot(a, fb, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "between variance",xlab = "Age", ylim=c(0,9),
		main = "",
		las = 1,
		cex.lab = 1.4)

text(c(15,20,15,20),fb[21,],c(1981,1991,2001,2011),pos=c(1,3,3,3)) 
dev.off()


# --------------------------------
# check Hal's version:
#.SD <- SCO[SCO$year == 1981 & SCO$sex == 1, ]
#Byrsex <- function(.SD,w="stationary"){
#	QXk <- acast(.SD, age ~quintile_2, value.var = "qx")
#	
#	# Hal has them stacked.
#	E <- c(apply(QXk, 2, getEta1k))
#	V <- c(apply(QXk, 2, getVk))
#	
#	pii <- apply(QXk, 2, function(qx){
#				Ui <- getNk(qx)
#				e1 <- rep(0,111)
#				e1[1] <- 1
#				ones <- rep(1,111)
#				t(ones)%*%(Ui %*% e1) / 5
#			})
#	
#	Ig <- diag(5)
#	eblank <- rep(0,111)
#	
#	for (i in 1:111){
#		epi    <- eblank
#		epi[i] <- 1
#		
#	}
#	
#	# left off here. odd notation.
#	
#	if (w == "raw"){
#		NXkts <- acast(.SD, age ~quintile_2, value.var = "N")
#		raww  <- NXkts / rowSums(NXkts)
#		out <- Vbetween(QXkts,raww)
#	}
#	out
#}
#Wyrsex <- function(.SD,w="stationary"){
#	QXkts <- acast(.SD, age ~quintile_2, value.var = "qx")
#	if (w == "stationary"){
#		out <- c(Vwithin(QXkts))
#	}
#	if (w == "raw"){
#		NXkts <- acast(.SD, age ~quintile_2, value.var = "N")
#		raww  <- NXkts / rowSums(NXkts)
#		out <- c(Vwithin(QXkts,raww))
#	}
#	out
#}
#
#










# end