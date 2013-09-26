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
	
	tb$par[T,A] <- TS
	tb$par[T,G] <- TW
	tb$par[C,G] <- TS
	tb$par[G,G] <- TW
	tb$par[G,T] <- TS
	tb$par[T,C] <- TW
	
	tb$apar[A,A] <- TS
	tb$apar[A,G] <- TW
	tb$apar[T,A] <- TS
	tb$apar[T,C] <- TW
	tb$apar[C,A] <- TW
	tb$apar[G,G] <- TS
	
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
	
	tb$par[T,A] <- IA
	tb$par[T,G] <- IA
	tb$par[C,G] <- IA
	tb$par[G,G] <- IB
	tb$par[G,T] <- IB
	tb$par[T,C] <- IB
	
	tb$apar[A,A] <- IC
	tb$apar[A,G] <- ID
	tb$apar[T,A] <- IC
	tb$apar[T,C] <- IE
	tb$apar[C,A] <- ID
	tb$apar[G,G] <- IE
	
	return(tb)
}
