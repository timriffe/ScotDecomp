
# Author: tim
###############################################################################

#SCO <- local(get(load("/home/tim/git/ScotDecomp/ScotDecomp/Data/SCOlong.Rdata")))
SCO <- read.csv("/home/tim/git/ScotDecomp/ScotDecomp/Data/lifetables_quintiles_scotland.csv")
SCO$Fac <- 3
SCO$Fac[SCO$year == 1991] <- 2

SCO$N <- SCO$N / SCO$Fac
SCO$D <- SCO$D / SCO$Fac

library(xtable)

# N per quintile, year, sex

Sums <- tapply(SCO$N, list(SCO$year, SCO$quintile_2, SCO$sex),sum)

# total pop
SumM <- t(Sums[,,1])
SumF <- t(Sums[,,2])

# e0
me0  <- acast(SCO[SCO$age==0 & SCO$sex == 1, ], quintile_2 ~ year, value.var = "ex")
fe0  <- acast(SCO[SCO$age==0 & SCO$sex == 2, ], quintile_2 ~ year, value.var = "ex")

# e35
me35 <- acast(SCO[SCO$age==35 & SCO$sex == 1, ], quintile_2 ~ year, value.var = "ex")
fe35 <- acast(SCO[SCO$age==35 & SCO$sex == 2, ], quintile_2 ~ year, value.var = "ex")

# sd0
# sd35
mid <- SCO$sex == 1
library(data.table)
SCO  <- data.table(SCO)
SCO  <- SCO[,list(V = c(getVk(qx)),age=unique(age)),by=list(year,quintile_2,sex)]
SCO  <- data.frame(SCO)

mV0  <- sqrt(acast(SCO[SCO$age==0 & SCO$sex == 1, ], quintile_2 ~ year, value.var = "V"))
fV0  <- sqrt(acast(SCO[SCO$age==0 & SCO$sex == 2, ], quintile_2 ~ year, value.var = "V"))

# e35
mV35 <- sqrt(acast(SCO[SCO$age==35 & SCO$sex == 1, ], quintile_2 ~ year, value.var = "V"))
fV35 <- sqrt(acast(SCO[SCO$age==35 & SCO$sex == 2, ], quintile_2 ~ year, value.var = "V"))

# --- now make interleaved dfs

# 12 rows:
MaleTab <- matrix(nrow=18,ncol=8)
MaleTab[seq(1,16,by=3),seq(1,7,2)]     <- round(SumM/1000)*1000
MaleTab[seq(1,16,by=3)+1,seq(1,7,2)]   <- round(me0,1)
MaleTab[seq(1,16,by=3)+2,seq(1,7,2)]   <- round(mV0,1)
MaleTab[seq(1,16,by=3)+1,seq(1,7,2)+1] <- round(me35,1)
MaleTab[seq(1,16,by=3)+2,seq(1,7,2)+1] <- round(mV35,1)

print(xtable(MaleTab,digits=1))
?xtable
# 5 dfs, 6 major rows,

#library(RColorBrewer)
#blues3 <- brewer.pal(7,"Blues")[-c(1,2)]
#red3   <- brewer.pal(7,"Reds")[-c(1,2)]




