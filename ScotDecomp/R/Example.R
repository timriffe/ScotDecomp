

# ok, decompose within vs between variance for 2 and N populations.

# 1) get some example data:
library(HMDHFDplus)
LT <- readHMDweb("SWE","mltper_1x1",us,pw)

head(LT)

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
	
	Eta1s <- apply(QXk, 2, getEta1k)
	(Eta1s ^ 2) %*% pik - (Eta1s %*% pik) ^ 2
}

getVk(qx)
plot(sqrt(getVk(qx)))

QXk <- cbind(qx, LT$qx[LT$Year == 2000])

Vwithin(QXk) +
Vbetween(QXk)

# how to average qx's together: use lx as weights.
lx1 <- getNk(QXk[,1])[, 1]
lx2 <- getNk(QXk[,2])[, 1]



cbind(Vwithin(QXk) - Vbetween(QXk),
getVk((lx1 * QXk[,1] + lx2 * QXk[,2] ) / (lx1 + lx2)),
getVk(rowMeans(QXk)))

plot(Vwithin(QXk) + Vbetween(QXk))
lines(getVk(rowMeans(QXk)))
lines(getVk((lx1 * QXk[,1] + lx2 * QXk[,2] ) / (lx1 + lx2)))
lines(Vwithin(QXk))
lines(Vbetween(QXk))