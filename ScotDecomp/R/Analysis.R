# Author: tim
###############################################################################
# TODO: condense into functions and iteration. Would be slicker code and easier to redo variants.

SCO <- read.csv("/home/tim/git/ScotDecomp/ScotDecomp/Data/lifetables_quintiles_scotland.csv")
source("/home/tim/git/ScotDecomp/ScotDecomp/R/Functions.R")
SCO <- SCO[SCO$quintile_2 != 999, ]
library(reshape2)
library(data.table)

# calculate things manually (too much typing sorry)

# make a function that calcs within and between for a year and sex:
# add raw weighting option as well
Byrsex <- function(.SD,w="stationary"){
	QXkts <- acast(.SD, age ~quintile_2, value.var = "qx")
	if (w == "stationary"){
		out <- Vbetween(QXkts)
	}
	out
}
Wyrsex <- function(.SD,w="stationary"){
	QXkts <- acast(.SD, age ~quintile_2, value.var = "qx")
	Vwithin(QXkts)
	if (w == "stationary"){
		out <- Vwithin(QXkts)
	}
	out
}
SCO <- data.table(SCO)
SCO[,B := Byrsex(.SD), by = list(year, sex)]
SCO[,W := Wyrsex(.SD), by = list(year, sex)]

# from here down not reworked yet
a <- 0:85

png("/home/tim/git/ScotDecomp/ScotDecomp/Figures/BetwweenAgeTime.png",width=800,height=400)
par(mfrow=c(1,2))
matplot(a, mp, type = 'l', col = gray(c(.6,.4,.2,0)),lwd = c(2,1.5,1.2,1),
		lty=1,
		ylab = "proportion",xlab = "age", ylim=c(0,.04),
		main = "Variance due to between-group differences, Males")
text(20,mp[21, ], c(1981,1991,2001,2011),pos=3,cex=.8)

abline(v=35
matplot(a, fp, type = 'l', col = gray(c(.6,.4,.2,0)),lwd = c(2,1.5,1.2,1),
		lty=1,
		ylab = "proportion",xlab = "age", ylim=c(0,.04),
		main = "Variance due to between-group differences, Females")
text(40,fp[41, ], c(1981,1991,2001,2011),pos=c(3,3,3,1),cex=.8)
abline(v=35)
dev.off()

# total variance:


# ---------------------------------------------
# this is a dumb repetition of the above using natural weights for pik
#

