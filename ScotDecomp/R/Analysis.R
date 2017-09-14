# Author: tim
###############################################################################


SCO <- read.csv("/home/tim/git/ScotDecomp/ScotDecomp/Data/lifetables_quintiles_scotland.csv")
source("/home/tim/git/ScotDecomp/ScotDecomp/R/Example.R")

library(reshape2)
QX1981m <- acast(SCO[SCO$year == 1981 & SCO$sex == 1 & SCO$quintile_2 != 999, ], age ~quintile_2, value.var = "qx")
QX1981f <- acast(SCO[SCO$year == 1981 & SCO$sex == 2 & SCO$quintile_2 != 999, ], age ~quintile_2, value.var = "qx")

QX1991m <- acast(SCO[SCO$year == 1991 & SCO$sex == 1 & SCO$quintile_2 != 999, ], age ~quintile_2, value.var = "qx")
QX1991f <- acast(SCO[SCO$year == 1991 & SCO$sex == 2 & SCO$quintile_2 != 999, ], age ~quintile_2, value.var = "qx")

QX2001m <- acast(SCO[SCO$year == 2001 & SCO$sex == 1 & SCO$quintile_2 != 999, ], age ~quintile_2, value.var = "qx")
QX2001f <- acast(SCO[SCO$year == 2001 & SCO$sex == 2 & SCO$quintile_2 != 999, ], age ~quintile_2, value.var = "qx")

QX2011m <- acast(SCO[SCO$year == 2011 & SCO$sex == 1 & SCO$quintile_2 != 999, ], age ~quintile_2, value.var = "qx")
QX2011f <- acast(SCO[SCO$year == 2011 & SCO$sex == 2 & SCO$quintile_2 != 999, ], age ~quintile_2, value.var = "qx")

# calculate things manually (too much typing sorry)

# 1981
B1981m <- Vbetween(QX1981m)
W1981m <- Vwithin(QX1981m)

B1981f <- Vbetween(QX1981f)
W1981m <- Vwithin(QX1981m)

# 1991
B1991m <- Vbetween(QX1991m)
W1991m <- Vwithin(QX1991m)

B1991f <- Vbetween(QX1991f)
W1991f <- Vwithin(QX1991f)

# 2001
B2001m <- Vbetween(QX2001m)
W2001m <- Vwithin(QX2001m)

B2001f <- Vbetween(QX2001f)
W2001f <- Vwithin(QX2001f)

# 2011
B2011m <- Vbetween(QX2011m)
W2011m <- Vwithin(QX2011m)

B2011f <- Vbetween(QX2011f)
W2011f <- Vwithin(QX2011f)

# prop between
mp <- cbind(B1981m / (B1981m + W1981m),
	  B1991m / (B1991m + W1991m),
	  B2001m / (B2001m + W2001m),
	  B2011m / (B2011m + W2011m)
)
fp <- cbind(B1981f / (B1981f + W1981f),
		B1991f / (B1991f + W1991f),
		B2001f / (B2001f + W2001f),
		B2011f / (B2011f + W2011f)
)

a <- 0:85

png("/home/tim/git/ScotDecomp/ScotDecomp/Figures/BetwweenAgeTime.png",width=800,height=400)
par(mfrow=c(1,2))
matplot(a, mp, type = 'l', col = gray(c(.6,.4,.2,0)),lwd = c(2,1.5,1.2,1),
		lty=1,
		ylab = "proportion",xlab = "age", ylim=c(0,.04),
		main = "Variance due to between-group differences, Males")
text(20,mp[21, ], c(1981,1991,2001,2011),pos=3,cex=.8)

matplot(a, fp, type = 'l', col = gray(c(.6,.4,.2,0)),lwd = c(2,1.5,1.2,1),
		lty=1,
		ylab = "proportion",xlab = "age", ylim=c(0,.04),
		main = "Variance due to between-group differences, Females")
text(40,fp[41, ], c(1981,1991,2001,2011),pos=c(3,3,3,1),cex=.8)
dev.off()

# total variance:

mt <- cbind((B1981m + W1981m),
		(B1991m + W1991m),
		(B2001m + W2001m),
		(B2011m + W2011m)
)
ft <- cbind((B1981f + W1981f),
		(B1991f + W1991f),
		(B2001f + W2001f),
		(B2011f + W2011f)
)


par(mfrow=c(1,2))
matplot(a, sqrt(mt), type = 'l', col = gray(c(.6,.4,.2,0)),lwd = c(2,1.5,1.2,1),
		lty=1,
		ylab = "sd",xlab = "age",ylim=c(0,16),
		main = "Total sd, Males")

matplot(a, sqrt(ft), type = 'l', col = gray(c(.6,.4,.2,0)),lwd = c(2,1.5,1.2,1),
		lty=1,
		ylab = "proportion",xlab = "age",ylim=c(0,16),
		main = "Total sd, Males")





