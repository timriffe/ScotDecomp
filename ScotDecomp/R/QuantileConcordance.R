
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
Car$quintile2011 <- 6 - Car$quintile2011
tail(Car)

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


tab1981 <- table(Car$quintile1981,mode)
tab1991 <- table(Car$quintile1991,mode)
tab2001 <- table(Car$quintile2001,mode)
tab2011 <- table(Car$quintile2011,mode)

prop_same <- function(tab){
	sum(diag(tab)) / sum(tab)
}
prop_same(tab1981) # 0.7534653
prop_same(tab1991) # 0.7712288
prop_same(tab2001) # 0.7930693
prop_same(tab2011) # 0.7440711

prop_off <- function(tab, offset = 0){
	ind <- abs(row(tab) - col(tab)) == offset
	sum(tab[ind]) / sum(tab)
}
prop_same(tab1981) + prop_off(tab1981,1) # 0.9257426
prop_same(tab1991) + prop_off(tab1991,1) # 0.9250749
prop_same(tab2001) + prop_off(tab2001,1) # 0.9693069
prop_same(tab2011) + prop_off(tab2011,1) # 0.958498





