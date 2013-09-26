###
## Triplex default score and isogroup tables
##
## Author: Jiri Hon
## Date: 2013/09/25
## Package: triplex
##


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
	
	tb$apar[A,T] <- TS
	tb$apar[A,C] <- TW
	tb$apar[T,T] <- TS
	tb$apar[T,G] <- TW
	tb$apar[C,T] <- TW
	tb$apar[G,C] <- TS
	
	return(tb)
}

###
## Get default isogroup table
##
triplex.group.table <- function()
{
	names <- c("A", "C", "G", "T")
	mx <- matrix(IN, NBASES, NBASES, dimnames=list(names, names))
	tb <- list(par=mx, apar=mx)
	
	tb$par[T,T] <- IA
	tb$par[T,C] <- IA
	tb$par[C,C] <- IA
	tb$par[G,C] <- IB
	tb$par[G,A] <- IB
	tb$par[T,G] <- IB
	
	tb$apar[A,T] <- IC
	tb$apar[A,C] <- ID
	tb$apar[T,T] <- IC
	tb$apar[T,G] <- IE
	tb$apar[C,T] <- ID
	tb$apar[G,C] <- IE
	
	return(tb)
}
