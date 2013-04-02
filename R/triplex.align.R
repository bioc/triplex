###
## Triplex align R interface
##
## Author: Jiri Hon
## Date: 2012/11/01
## Package: triplex
##


###
## Get DNA complementary strand
##
complement <- function(seq)
{
	chartr("atgcATGC", "tacgTACG", seq)
}


###
## Align triplex
##
triplex.align <- function(t)
{
	if (class(t) != "TriplexViews")
		stop("Triplex must be TriplexViews object")
	
	if (length(t) != 1)
		stop("TriplexViews object must be of length 1")
	
	seq <- subseq(subject(t), start(t), end(t))
	
	p <- t@params
	p[MIN_LOOP] <- to_double(lwidth(t))
	#p[MIN_LOOP] <- p[MIN_LOOP] + 1
	type <- to_double(type(t))
	
	triplex <- .Call("triplex_align", seq, type, p)
	
	# Make reverse complement
	if (strand(t) == "-")
		triplex <- reverse(complement(triplex))
	
	return(triplex)
}


###
## Convert triplex alignment into DNAStringSet
##
alignment.DNAStringSet <- function(alignment, type)
{
	n <- nchar(alignment)
	seq <- substring(alignment, 1:n, 1:n)
	loop <- which(seq == "=")
	
	seq.left <- paste(seq[1:(loop[1] - 1)], collapse='')
	seq.loop <- paste(seq[(loop[1] + 1):(loop[2] - 1)], collapse='')
	seq.right <- paste(seq[(loop[2] + 1):n], collapse='')
	
	if (type == 0)
	{
		a <- DNAStringSet(c(
			seq.left,
			complement(seq.left),
			reverse(seq.right),
			seq.loop))
		names(a) <- c(
			"plus",
			"minus",
			"para-plus",
			"loop")
	}
	else if (type == 1)
	{
		a <- DNAStringSet(c(
			seq.left,
			reverse(complement(seq.right)),
			reverse(seq.right),
			seq.loop))
		names(a) <- c(
			"para-plus",
			"minus",
			"plus",
			"loop")
	}
	else if (type == 2)
	{
		a <- DNAStringSet(c(
			seq.left,
			reverse(complement(seq.right)),
			reverse(seq.right),
			seq.loop))
		names(a) <- c(
			"para-minus",
			"plus",
			"minus",
			"loop")
	}
	else if (type == 3)
	{
		a <- DNAStringSet(c(
			seq.left,
			complement(seq.left),
			reverse(seq.right),
			seq.loop))
		names(a) <- c(
			"minus",
			"plus",
			"para-minus",
			"loop")
	}
	else if (type == 4)
	{
		a <- DNAStringSet(c(
			reverse(complement(seq.right)),
			reverse(seq.right),
			seq.left,
			seq.loop))
		names(a) <- c(
			"plus",
			"minus",
			"anti-minus",
			"loop")
	}
	else if (type == 5)
	{
		a <- DNAStringSet(c(
			reverse(seq.right),
			seq.left,
			complement(seq.left),
			seq.loop))
		names(a) <- c(
			"anti-minus",
			"minus",
			"plus",
			"loop")
	}
	else if (type == 6)
	{
		a <- DNAStringSet(c(
			reverse(seq.right),
			seq.left,
			complement(seq.left),
			seq.loop))
		names(a) <- c(
			"anti-plus",
			"plus",
			"minus",
			"loop")
	}
	else
	{
		a <- DNAStringSet(c(
			reverse(complement(seq.right)),
			reverse(seq.right),
			seq.left,
			seq.loop))
		names(a) <- c(
			"minus",
			"plus",
			"anti-plus",
			"loop")
	}
	return(a)
}


###
## Get triplex alignment for further processing
##
triplex.alignment <- function(triplex)
{
	alignment <- triplex.align(triplex)
	type <- type(triplex)
	alignment.DNAStringSet(alignment, type)
}
