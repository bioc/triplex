###
## Triplex align R interface
##
## Author: Jiri Hon
## Date: 2013/09/25
## Package: triplex
##

NBASES = 4

TM = -9
TS = 2
TW = 1

A = 1
C = 2
G = 3
T = 4


###
## Get default score tables
##
triplex.score.table <- function()
{
	names <- c("A", "C", "G", "T")
	mx <- matrix(TM, NBASES, NBASES, dimnames=list(names, names))
	tb <- list(par=mx, apar=mx)
	
	tb$par[T,T] <- TS
	tb$par[T,C] <- TW
	tb$par[C,C] <- TS
	tb$par[G,C] <- TW
	tb$par[G,A] <- TS
	tb$par[T,G] <- TW
	
	return(tb)
}
