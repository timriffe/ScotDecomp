
# ---------------------
source(file.path("R","0_Functions.R"))
#source("R/Extrapolate.R")
# make a function that calcs within and between for a year and sex:
# raw weighting not prepared because N in orig dat aonly up to 85+,
# could resistribute 85+ according to l(85+), but that's for another
# day. 
Byr <- function(.SD,w="stationary"){
	# only difference is that we have sex in columns...
	QXkts <- acast(.SD, age ~ sex, value.var = "qx")
	if (w == "stationary"){
		out <- Vbetween(QXkts)
	}
	if (w == "raw"){
		# only difference is that we have sex in columns...
		NXkts <- acast(.SD, age ~ sex, value.var = "N")
		raww  <- NXkts / rowSums(NXkts)
		out   <- Vbetween(QXkts, raww)
	}
	out
}
Wyr <- function(.SD,w="stationary"){
	# only difference is that we have sex in columns...
	QXkts <- acast(.SD, age ~ sex, value.var = "qx")
	if (w == "stationary"){
		out <- c(Vwithin(QXkts))
	}
	if (w == "raw"){
		# only difference is that we have sex in columns...
		NXkts <- acast(.SD, age ~ sex, value.var = "N")
		raww  <- NXkts / rowSums(NXkts)
		out   <- c(Vwithin(QXkts, raww))
	}
	out
}

SCO <- local(get(load(file.path("Data","Derived","SCOlong.Rdata"))))
# just total pop, so we have a strata for sex only due to this selection.
SCO <- SCO[SCO$quintile_2 == 999, ]
# calculate things manually (too much typing sorry)
SCO <- data.table(SCO)

# when we did the decomp by quintiles, the by = list()
# argument included sex and year, but since it's over sex
# and there are no quintiles then it's over year only.
SCOB <- SCO[,list(age = unique(age),
				Bst = Byr(.SD), 
				Wst = Wyr(.SD)), 
		by = list(year)]

# derive total, proportion between, and standard deviation
SCOB$V     <- SCOB$B + SCOB$W
SCOB$propB <- SCOB$B / SCOB$V
SCOB$sd    <- sqrt(SCOB$V)

# save out results to csv
write.csv(SCOB,file=file.path("Data","Derived","SexDecomp.csv"))


do.this <- FALSE
if (do.this){
#  a couple optional figures
probB <- acast(SCOB, age~year, value.var = "propB")

png("Figures/BetweenSexDecrease.png")
matplot(0:110,probB,type = 'l',col=gray(seq(0,.6,length=4)),
		main = "Proportion of total variance due to between-sex differences",xlab="Age",ylab = "Proportion")
text(20,probB["20", ], colnames(probB),pos=1)
dev.off()


absB <- acast(SCOB, age~year, value.var = "Bst")
png("Figures/BetweenSexAbs.png")
matplot(0:110,absB,type = 'l',col=gray(seq(0,.6,length=4)),
		main = "Total variance due to between-sex differences",xlab="Age",ylab = "Variance Component")
text(20,absB["20", ], colnames(absB),pos=1)
dev.off()
}
#absW <- acast(SCOB, age~year, value.var = "Wst")
##png("Figures/BetweenSexAbs.png")
#matplot(0:110,absW,type = 'l',col=gray(seq(0,.6,length=4)),
#		main = "Proportion of total variance due to between-sex differences",xlab="Age",ylab = "Proportion")
#text(20,absW["20", ], colnames(absW),pos=1)
#dev.off()


# another comparison: what do we get if we decompose wrt ONLY quintiles 1 vs 5?
# i.e. throw out 2,3,4.
# ---------------------------------------
# same functions as in Analysis
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

SCO <- local(get(load(file.path("Data","Derived","SCOlong.Rdata"))))
# select down
SCO <- SCO[SCO$quintile_2 %in% c(1,5), ]

SCO <- data.table(SCO)

# this might look familiar by now
SCOB <- SCO[,list(age = unique(age),
				Bst = Byrsex(.SD), 
				Wst = Wyrsex(.SD)), 
		by = list(year,sex)]


SCOB$V     <- SCOB$B + SCOB$W
SCOB$propB <- SCOB$B / SCOB$V
SCOB$sd    <- sqrt(SCOB$V)

# save out results to csv
write.csv(SCOB,file=file.path("Data","Derived","ExtremeDecomp.csv"))


if (do.this){
	# get some matrices to plot
	probBm <- acast(SCOB[SCOB$sex == 1, ], age~year, value.var = "propB")
	probBf <- acast(SCOB[SCOB$sex == 2, ], age~year, value.var = "propB")
# as high as 8%!
png("Figures/BetweenPropExtremeMales.png")
matplot(0:110,probBm,type = 'l',col=gray(seq(0,.6,length=4)),
		main = "Proportion of subpopulation [1,5] variance due to group differences, Males",xlab="Age",ylab = "Proportion")
text(20,probBm["20", ], colnames(probBm),pos=1)
dev.off()

# as high as 4%
png("Figures/BetweenPropExtremeFemales.png")
matplot(0:110,probBf,type = 'l',col=gray(seq(0,.6,length=4)),
		main = "Proportion of subpopulation [1,5] variance due to group differences, Females",xlab="Age",ylab = "Proportion")
text(20,probBf["20", ], colnames(probBf),pos=1)
dev.off()

}