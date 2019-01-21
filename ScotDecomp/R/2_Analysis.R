
# ---------------------
source(file.path("R","0_Functions.R"))

SCO <- local(get(load(file.path("Data","Derived","SCOlong.Rdata"))))

# remove the Scotland aggregate. It isn't the same as
# the synthetic blend of quintile lifetables. cuz the
# weights would be different.
SCO <- SCO[SCO$quintile_2 != 999, ]

# convert data type for easier calcs
SCO <- data.table(SCO)

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
save(SCOB,file=file.path("Data","Derived","SCOB.Rdata"))



# Figure production in R (not used in manuscript)

do.this <- FALSE
if (do.this){
mp         <- acast(SCOB[SCOB$sex == 1, ], age~year, value.var = "propB")
fp         <- acast(SCOB[SCOB$sex == 2, ], age~year, value.var = "propB")
mp         <- mp[1:86, ]
fp         <- fp[1:86, ]
# from here down not reworked yet
a          <- 0:85

maxA <- function(x,age=1:length(x)-1,trunc = 80){
	which.max(x[age<trunc])-1
}
malemax   <- tapply(SCOB$propB[SCOB$sex == 1 ],SCOB$year[SCOB$sex == 1 ],maxA)
femalemax <- tapply(SCOB$propB[SCOB$sex == 2 ],SCOB$year[SCOB$sex == 2 ],maxA)


graphics.off()
pdf("Figures/BetweenPropMales.pdf")
matplot(a, mp, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "Proportion",xlab = "Age", ylim=c(0,.04),
		main = "",
		las = 1,
		cex.lab = 1.4)
#abline(v=35)
text(10,mp[11, ], c(1981,1991,2001,2011),pos=3,cex=1)
for (i in 1:4){
	points(malemax[i],mp[malemax[i]+1, i], pch = 16)
	text(malemax[i],mp[malemax[i]+1, i],malemax[i],pos=1)
}

dev.off()

pdf("Figures/BetweenPropMalesRelChg.pdf")
percChgm <- 100*t(mp[c(1,31), ] / mp[c(1,31), 1])
matplot(c(1981,1991,2001,2011),percChgm, type = 'o', 
		ylab = "relative change (%) in prop between",
		xlab = "Year",
		ylim=c(100,280),
		las = 1,
		pch= 16)
text(2001, percChgm[3,]-20, c("Age 0", "Age 30"), pos = 1)
dev.off()

pdf("Figures/BetweenPropFemalesRelChg.pdf")
percChgf <- 100*t(fp[c(1,31), ] / fp[c(1,31), 1])
matplot(c(1981,1991,2001,2011),percChgf, type = 'o', 
		ylab = "relative change (%) in prop between",
		xlab = "Year",
		ylim=c(100,280),
		las = 1,
		pch= 16)
text(2001, percChgf[3,]-5, c("Age 0", "Age 30"), pos = 1)
dev.off()

pdf("Figures/BetweenPropFemales.pdf")
matplot(a, fp, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "Proportion",xlab = "Age", ylim=c(0,.04),
		main = "",
		las = 1,
		cex.lab = 1.4)
text(10,fp[11, ], c(1981,1991,2001,2011),pos=c(1,3,3,3),cex=1)
for (i in 1:4){
	points(femalemax[i],fp[femalemax[i]+1, i], pch = 16)
	text(femalemax[i],fp[femalemax[i]+1, i],femalemax[i],pos=1)
}
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

# new graphs as of 23-Jan-2018
# absolute var within:
mw         <- acast(SCOB[SCOB$sex == 1, ], age~year, value.var = "Wst")
fw         <- acast(SCOB[SCOB$sex == 2, ], age~year, value.var = "Wst")
mw         <- mw[1:86, ]
fw         <- fw[1:86, ]
mv         <- acast(SCOB[SCOB$sex == 1, ], age~year, value.var = "V")
fv         <- acast(SCOB[SCOB$sex == 2, ], age~year, value.var = "V")
mv         <- mv[1:86, ]
fv         <- fv[1:86, ]
pdf("Figures/TotalwstMales.pdf")
matplot(a, mw, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "within variance",xlab = "Age", #ylim = c(0,9),
		main = "",
		las = 1,
		cex.lab = 1.4,
		ylim = c(0,260))
text(20,mw[21,],c(1981,1991,2001,2011)) 
dev.off()

pdf("Figures/TotalwstFemales.pdf")
matplot(a, fw, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "within variance",xlab = "Age", #ylim=c(0,9),
		main = "",
		las = 1,
		cex.lab = 1.4,
		ylim = c(0,260))

text(c(15,20,15,20),fw[21,],c(1981,1991,2001,2011),pos=c(1,3,3,3)) 
dev.off()

# proportion within
pdf("Figures/WithinPropMales.pdf")
matplot(a, mw/mv, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "Proportion",xlab = "Age", ylim = c(.96,1),
		main = "",
		las = 1,
		cex.lab = 1.4)
		#ylim = c(0,260))
text(10,mw[11,]/mv[11,],c(1981,1991,2001,2011),pos=3) 
dev.off()

pdf("Figures/WithinPropFemales.pdf")
matplot(a, fw/fv, type = 'l', col = gray(c(.7,.5,.3,0)),lwd = c(3,2,1.5,1),
		lty=1,
		ylab = "Proportion",xlab = "Age", ylim=c(.96,1),
		main = "",
		las = 1,
		cex.lab = 1.4)
		#ylim = c(0,260))

text(10,fw[11,]/fv[11,],c(1981,1991,2001,2011),pos=c(3,1,1,1)) 
dev.off()

}



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