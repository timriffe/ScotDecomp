# this is a data check, independent of th rest of analyses

source(file.path("R","0_Functions.R"))
# data path valid only if on MPIDR PC:


Car <- read.dta(file.path("Data","Inputs","Carstairs_V12.dta"))
Car <- data.frame(Car)

# flip 2011 quintiles due to odd coding switch
Car$quintile[Car$year == 2011] <- 6 - Car$quintile[Car$year == 2011]

# how many unique codes per year
length(unique(Car$pcsector[Car$year == 1981]))
length(unique(Car$pcsector[Car$year == 1991]))
length(unique(Car$pcsector[Car$year == 2001]))
length(unique(Car$pcsector[Car$year == 2011]))

# how many cases of unsplit vs twice split bla bla per year
tab1981 <- table(Car$pcsector[Car$year == 1981])
tab1991 <- table(Car$pcsector[Car$year == 1991])
tab2001 <- table(Car$pcsector[Car$year == 2001])
tab2011 <- table(Car$pcsector[Car$year == 2011])

table(tab1981)
table(tab1991)
table(tab2001)
table(tab2011)

# -------------------------------
# codes of postal code sectors that never switch
keep1981 <- names(tab1981)[tab1981 == 1]
keep1991 <- names(tab1991)[tab1991 == 1]
keep2001 <- names(tab2001)[tab2001 == 1]
keep2011 <- names(tab2011)[tab2011 == 1]

Car$once <- FALSE
Car$once[Car$year == 1981 & Car$pcsector %in% keep1981] <- TRUE
Car$once[Car$year == 1991 & Car$pcsector %in% keep1991] <- TRUE
Car$once[Car$year == 2001 & Car$pcsector %in% keep2001] <- TRUE
Car$once[Car$year == 2011 & Car$pcsector %in% keep2011] <- TRUE

# ca 180 / census are split
#sum(!Car$once)/4
# --------------------------------


Quint1 <- acast(Car[Car$once, ], pcsector~year, value.var = "quintile")
mode   <- apply(Quint1,1,get_modal)

(tab1981 <- table(Quint1[,1],mode))
(tab1991 <- table(Quint1[,2],mode))
(tab2001 <- table(Quint1[,3],mode))
(tab2011 <- table(Quint1[,4],mode))

prop_same(tab1981) # 0.7580438
prop_same(tab1991) # 0.7959698
prop_same(tab2001) # 0.8085352
prop_same(tab2011) # 0.752

prop_same(tab1981) + prop_off(tab1981,1) # 0.9420849
prop_same(tab1991) + prop_off(tab1991,1) # 0.9571788
prop_same(tab2001) + prop_off(tab2001,1) # 0.9861592
prop_same(tab2011) + prop_off(tab2011,1) # 0.9714286

# --------------------------------------
do.this <- FALSE
if (do.this){
# still not convinced? How about we look at deprivation score
# distributions: observed quintile vs modal quintile, per year


z1 <- acast(Car[Car$once, ], pcsector~year, value.var = "carstair")
breaks <- seq(-10,10,by=.2)

compare <- function(z1,quint,mode,col1 = "#0000FF50", col2 = "#FF000050", breaks = seq(-10,10,by=.2),year=1,q=1){
	hist(z1[quint[,year] == q,year],breaks=breaks, col = col1)
	hist(z1[mode == q,year], add = TRUE, col = col2,breaks=breaks)
}
breaks <- seq(-16,16,by=.5)
par(mfrow=c(4,5), mai=c(.1,.1,.1,.1))
for (yr in 1:4){
     for (quint in 1:5){
			compare(z1,Quint1,mode,breaks=breaks,year=yr,q=quint)
	}
}
}





