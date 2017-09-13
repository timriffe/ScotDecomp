

# ok, decompose within vs between variance for 2 and N populations.

# 1) get some example data:
library(HMDHFDplus)
LT <- readHMDweb("SWE","mltper_1x1",us,pw)


qx <- LT$qx[1:111]

# first eq in Hal's notes,
# we want px aka (1-qx) in subdaig
getUk <- function(qx){
	N <- length(qx)
	# diag
	O <- diag(1 - qx[-N])
	# make subdiag
	U <- cbind(rbind(0, O), 0)
	# closeout
	U[N, N] <- 1 - qx[N]
	U
}

# get N from qx
getNk <- function(qx){
	Uk <- getUk(qx)
	Nk <- solve(diag(length(qx))-Uk)
	Nk
}

# get eta 1 from qx
# eta 1 is close to ex.
# the difference between eta 1 and lifetable ex is that eta1
# assumes survival until the end of the age interval, so it's
# a bit higher.
getEta1k <- function(qx){
	N    <- length(qx)
	Nk   <- getNk(qx)
	ones <- rep(1,N)
	t(t(ones) %*% Nk)
}
# get eta 2 from qx
getEta2k <- function(qx){
	N    <- length(qx)
	Nk   <- getNk(qx)
	ones <- rep(1,N)
	t(t(ones) %*% Nk %*% (2 * Nk - diag(N)))
}

# the variance of remaining lifespan
getVk  <- function(qx){
	Eta1 <- getEta1k(qx)
	Eta2 <- getEta2k(qx)
	Eta2 - Eta1 ^ 2
}

# calculate within-variance (QXk is age in rows and pops in columns)
# length(pik) must equal ncol(QXk)
Vwithin <- function(QXk, pik){
	# if pi not given, let's assume it's
	# uniform
	if (missing(pik)){
		pik <- rep(1 / ncol(QXk), ncol(QXk))
	}
	stopifnot(length(pik) == ncol (QXk))
	
	VX <- apply(QXk, 2, getVk)
	
	Vw <- VX %*% pik
	Vw
}

# calculate between-variance(QXk is age in rows and pops in columns)
# length(pik) must equal ncol(QXk)
Vbetween <- function(QXk, pik){
	# if pi not given, let's assume it's
	# uniform
	if (missing(pik)){
		pik <- rep(1 / ncol(QXk), ncol(QXk))
	}
	stopifnot(length(pik) == ncol (QXk))
	# make sure sums to 1
	pik <- pik / sum(pik)
	
	Eta1s <- apply(QXk, 2, getEta1k)
	((Eta1s ^ 2) %*% pik) - ((Eta1s %*% pik) ^ 2)
}

# tests: what do within and between variance sum to? total?
QXk <- cbind(qx, LT$qx[LT$Year == 2000])

Vwithin(QXk) +
Vbetween(QXk)
Vbetween(QXk) + Vwithin(QXk)



library(DistributionTTD)
qx2dx <- function(qx){
	N <- length(qx)
	lx <- cumprod(c(1,1-qx))
	dx <- -diff(c(lx,0))
	dx[N] <- lx[N]
	dx[1:N]
}
dx <- qx2dx(qx)
# holy shit they're identical
cbind(getVk(qx),
momentN(dx, n=2))

# so hwo to combine pops:
dx1   <- qx2dx(QXk[,1])
dx2   <- qx2dx(QXk[,2])
dxtot <- dx1 + dx2
dxtot <- dxtot/sum(dxtot)

plot(dx1,col = "red",type='l')
lines(dx2,col = "blue")
lines(dxtot,col="purple")


v1   <- getVk(dx2qx(dx1))
v2   <- getVk(dx2qx(dx2))
vtot <- getVk(dx2qx(dxtot))

plot(vtot,col="purple",type='l')
lines(v1,col = "red")
lines(v2,col = "blue")
lines(Vbetween(QXk) + Vwithin(QXk))
v1_2   <- momentN(dx1,2)
v2_2   <- momentN(dx2,2)
vtot_2 <- momentN(dxtot,2)
plot(vtot_2,col="purple",type='l')
lines(v1_2,col = "red")
lines(v2_2,col = "blue")


plot(momentN(dxtot, n=2))
lines(Vbetween(QXk) + Vwithin(QXk))

# dx2qx:
dx2qx<- function(dx){
	lx <- rev(cumsum(rev(dx)))
	dx / lx
}

plot(Vbetween(QXk) + Vwithin(QXk))

# 3 identical quantities:
lines(getVk(dx2qx(dxtot)))
lines(momentN(dxtot, n=2),lty=2,col="red",lwd=2)
lx1 <- qx2lx(QXk[,1])/1e5;lx1 <- qx2lx(QXk[,2])/1e5
lines(getVk((lx1 * QXk[,1] + lx2 * QXk[,2])/(lx1+lx2)),col = "blue",lty=3)
