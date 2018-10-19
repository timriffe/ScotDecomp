
me <- system("whoami",intern=TRUE)
if (me == "tim"){
	setwd("/home/tim/git/ScotDecomp/ScotDecomp")
} 
if (me == "mpidr_d\\riffe"){
	setwd("U:/git/ScotDecomp/ScotDecomp")
}

library("foreign")

# data path valid only if on MPIDR PC:

data.path <- "N:/Rosie/For Tim/BMJOpen/Data/V12"

Car <- read.dta(file.path(data.path,"Carstairs_V12.dta"))
Car <- data.frame(Car)


get_modal <- function(x){
	tab <- table(unlist(x))
	as.integer(names(tab)[which.max(tab)])
}
mode <- apply(Car,1,get_modal)

Modes <- data.frame(pcsector = Car$pcsector, modal_quintile = mode)
write.csv(Modes, file = file.path(data.path,"CarMode.csv"))

# once-off spot check of heavy switchers
carvar <- rowSums((Car[,-1] - rowMeans(Car[,-1],na.rm=TRUE))^2,na.rm=TRUE) / rowSums(!is.na(Car[,-1]))
Car$carvar <- carvar
checkind <- rowSums(!is.na(Car[,2:5])) == 4
carcheck <- Car[checkind, ]
tail(carcheck[order(carcheck$carvar), ], 10)[,1]
# give Rosie the top 10 variance pcsectors, emailed




