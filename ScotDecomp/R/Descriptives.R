
# Author: tim
###############################################################################
# header to determine working directory
me <- system("whoami")
if (me == "tim"){
	setwd("/home/tim/git/ScotDecomp/ScotDecomp")
}
if (me == "mpidr_d\\seaman"){
	setwd("U:\\Conferences\\PAA\\2018 Denver\\within and between\\ScotDecomp")
}

# ---------------------
library(data.table)
library(xtable)
source("R/Functions.R")


#SCO <- local(get(load("/home/tim/git/ScotDecomp/ScotDecomp/Data/SCOlong.Rdata")))
SCO     <- read.csv("Data/lifetables_quintiles_scotland.csv")
SCOlong <- local(get(load("Data/SCOlong.Rdata")))
head(SCO)

# unadjust population size
SCO$Fac <- 3
SCO$Fac[SCO$year == 1991] <- 2

SCO$N   <- SCO$N / SCO$Fac
SCO$D   <- SCO$D / SCO$Fac



SCOlong <- data.table(SCOlong)
SCOlong <- SCOlong[,list(V = c(getVk(qx)), age = unique(age)), by = list(year,sex,quintile_2)]
SCOlong <- data.frame(SCOlong)

SCOlong <- SCOlong[SCOlong$age < 86, ]
SCO     <- merge(SCO, SCOlong)
SCO$sd  <- sqrt(SCO$V)

SCOcsv  <- SCO[SCO$age %in% c(0,35), ]
SCOcsv  <- SCOcsv[with(SCOcsv, order(year, sex, quintile_2,age)), ]

# save to csv
write.csv(SCOcsv, file = "Data/descriptives.csv", row.names = FALSE)
# -----------------------------------
# pop totals: just say in text ca 500k exposure for each individual lifetable.


# N per quintile, year, sex
#
#Sums <- tapply(SCO$N, list(SCO$year, SCO$quintile_2, SCO$sex),sum)
#
## total pop
#SumM <- t(Sums[,,1])
#SumF <- t(Sums[,,2])
#
## e0
#me0  <- acast(SCO[SCO$age==0 & SCO$sex == 1, ], quintile_2 ~ year, value.var = "ex")
#fe0  <- acast(SCO[SCO$age==0 & SCO$sex == 2, ], quintile_2 ~ year, value.var = "ex")
#
## e35
#me35 <- acast(SCO[SCO$age==35 & SCO$sex == 1, ], quintile_2 ~ year, value.var = "ex")
#fe35 <- acast(SCO[SCO$age==35 & SCO$sex == 2, ], quintile_2 ~ year, value.var = "ex")
#
## sd0
## sd35
#mid <- SCO$sex == 1
#library(data.table)
#SCO  <- data.table(SCO)
#SCO  <- SCO[,list(V = c(getVk(qx)),age=unique(age)),by=list(year,quintile_2,sex)]
#SCO  <- data.frame(SCO)
#
#mV0  <- sqrt(acast(SCO[SCO$age==0 & SCO$sex == 1, ], quintile_2 ~ year, value.var = "V"))
#fV0  <- sqrt(acast(SCO[SCO$age==0 & SCO$sex == 2, ], quintile_2 ~ year, value.var = "V"))
#
## e35
#mV35 <- sqrt(acast(SCO[SCO$age==35 & SCO$sex == 1, ], quintile_2 ~ year, value.var = "V"))
#fV35 <- sqrt(acast(SCO[SCO$age==35 & SCO$sex == 2, ], quintile_2 ~ year, value.var = "V"))
#
## --- now make interleaved dfs
#
## 12 rows:
##MaleTab <- matrix(nrow=18,ncol=8)
##MaleTab[seq(1,16,by=3),seq(1,7,2)]     <- round(SumM/1000)*1000
##MaleTab[seq(1,16,by=3)+1,seq(1,7,2)]   <- round(me0,1)
##MaleTab[seq(1,16,by=3)+2,seq(1,7,2)]   <- round(mV0,1)
##MaleTab[seq(1,16,by=3)+1,seq(1,7,2)+1] <- round(me35,1)
##MaleTab[seq(1,16,by=3)+2,seq(1,7,2)+1] <- round(mV35,1)
#
#
##print(xtable(cbind(round(SumM/1000,digits=0)),digits=0,caption="Population Size"))
##print(xtable(cbind(me0,fe0),digits=1,caption="Life Expectancy at birth"))
##print(xtable(cbind(me35,fe35),digits=1,caption="Life Expectancy at age 35"))
##print(xtable(cbind(mV0,fV0),digits=1,caption="Standard deviation, age 0"))
##print(xtable(cbind(mV35,mV35),digits=1,caption="Standard deviation, age 35"))
#
#
## 5 dfs, 6 major rows,
#
##library(RColorBrewer)
##blues3 <- brewer.pal(7,"Blues")[-c(1,2)]
##red3   <- brewer.pal(7,"Reds")[-c(1,2)]
#
#
#
#
