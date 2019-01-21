
# Author: tim
###############################################################################


SCO <- read.csv("/home/tim/git/ScotDecomp/ScotDecomp/Data/lifetables_quintiles_scotland.csv")

struct <- reshape2::acast(SCO[SCO$year == 2011 & SCO$quintile_2 != 999, ], age~quintile_2, sum, value.var = "N")

cols <- RColorBrewer::brewer.pal(7,"PuRd")[-c(1:2)]
png("ScotlandDep.png")
matplot(0:84,struct[-86,],
		type = 'l', col = cols,
		lwd = seq(3,1,length=5),
		lty = 1,
		las = 1,
		ylab = "Population",
		xlab = "Age",
		main = "Age by deprivation quintile Scotland, 2011",
		sub = "HT Rosie Seaman")
text(24,struct[25,],1:5,pos=3)
dev.off()
getwd()