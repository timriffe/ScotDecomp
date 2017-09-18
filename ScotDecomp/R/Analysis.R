# Author: tim
###############################################################################
# TODO: condense into functions and iteration. Would be slicker code and easier to redo variants.


source("/home/tim/git/ScotDecomp/ScotDecomp/R/Functions.R")
source("/home/tim/git/ScotDecomp/ScotDecomp/R/Extrapolate.R")

SCO <- local(get(load("/home/tim/git/ScotDecomp/ScotDecomp/Data/SCOlong.Rdata")))
SCO <- SCO[SCO$quintile_2 != 999, ]

library(reshape2)
library(data.table)

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

mp         <- acast(SCOB[SCOB$sex == 1, ], age~year, value.var = "propB")
fp         <- acast(SCOB[SCOB$sex == 2, ], age~year, value.var = "propB")
mp         <- mp[1:86, ]
fp         <- fp[1:86, ]
# from here down not reworked yet
a          <- 0:85

graphics.off()
pdf("/home/tim/git/ScotDecomp/ScotDecomp/Figures/BetweenPropMales.pdf")
matplot(a, mp, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "Proportion",xlab = "Age", ylim=c(0,.04),
		main = "",
		las = 1,
		cex.lab = 1.4)
abline(v=35)
text(20,mp[21, ], c(1981,1991,2001,2011),pos=3,cex=1)
dev.off()

pdf("/home/tim/git/ScotDecomp/ScotDecomp/Figures/BetweenPropFemales.pdf")
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

pdf("/home/tim/git/ScotDecomp/ScotDecomp/Figures/TotalsdMales.pdf")
matplot(a, msd, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "sd",xlab = "Age", ylim = c(0,16),
		main = "",
		las = 1,
		cex.lab = 1.4)
text(20,msd[21, "1981" ],1981,pos=2) 
text(30,msd[31, "2011" ],2011,pos=4) 
dev.off()

pdf("/home/tim/git/ScotDecomp/ScotDecomp/Figures/TotalsdFemales.pdf")
matplot(a, fsd, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "sd",xlab = "Age", ylim=c(0,16),
		main = "",
		las = 1,
		cex.lab = 1.4)

text(70,fsd[71, "1981"], 1981, pos=2)
text(70,fsd[71, "2011"], 2011, pos=4)
dev.off()

# end