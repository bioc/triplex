###
## TriplexViews class
##
## The TriplexViews class is the container for storing a set of triplex views
## (start/end locations, score, pvalue, etc.) on the same DNAString object,
## called the "subject" string.
##
## Author: Jiri Hon
## Date: 2012/10/21
## Package: triplex
##

setClass("TriplexViews",
	contains = "XStringViews",
	representation(
		params="numeric",
		score_table="list",
		group_table="list"
	)
)


###
## User-friendly constructor
##
TriplexViews <- function(
	subject, start, end, score, pvalue, insdel,
	type, lstart, lend, strand, params,
	score_table, group_table)
{
	view <- newViews(
		subject,
		start = start, end = end,
		Class = "TriplexViews"
	)
	view@params <- params
	view@score_table <- score_table
	view@group_table <- group_table
	
	mcols(view) <- DataFrame(
		score  = score,
		pvalue = pvalue,
		insdel = insdel,
		type   = type,
		lstart = lstart,
		lend   = lend,
		strand = strand
	)
	return(view)
}


###
## Set accessors to element metadata
##
setMethod("score", "TriplexViews", function(x) mcols(x)$score)

setGeneric("pvalue", function(x, ...) standardGeneric("pvalue"))
setMethod("pvalue", "TriplexViews", function(x) mcols(x)$pvalue)

setGeneric("ins", function(x, ...) standardGeneric("ins"))
setMethod("ins", "TriplexViews", function(x) mcols(x)$ins)

setMethod("type", "TriplexViews", function(x) mcols(x)$type)

setGeneric("lstart", function(x, ...) standardGeneric("lstart"))
setMethod("lstart", "TriplexViews", function(x) mcols(x)$lstart)

setGeneric("lend", function(x, ...) standardGeneric("lend"))
setMethod("lend", "TriplexViews", function(x) mcols(x)$lend)

setMethod("strand", "TriplexViews", function(x) mcols(x)$strand)

setGeneric("lwidth", function(x, ...) standardGeneric("lwidth"))
setMethod("lwidth", "TriplexViews", function(x)
{
	mcols(x)$lend - mcols(x)$lstart + 1
})


###
## The 2 helper functions below convert a given view on an XString object
## into a character-string.
##
## Both assume that 'start' <= 'end' (so they don't check it) and
## padd the result with spaces to produce the "margin effect"
## if 'start' or 'end' are out of limits.
##
## NOTE: Copied from Biostrings package, file XStringViews-class.R
##

###
## nchar(get_view(x, start, end)) is always end-start+1
##
get_view <- function(x, start, end, strand)
{
	s <- DNAString(subseq(x, start, end))
	if (strand == "-") s <- reverseComplement(s)
	as.character(s)
}

###
## nchar(get_snippet(x, start, end, snippetWidth)) is <= snippetWidth
##
get_snippet <- function(x, start, end, snippetWidth, strand)
{
	if (snippetWidth < 7)
		snippetWidth <- 7
	width <- end - start + 1
	if (width <= snippetWidth) {
		get_view(x, start, end, strand)
	} else {
		w1 <- (snippetWidth - 2) %/% 2
		w2 <- (snippetWidth - 3) %/% 2
		paste(get_view(x, start, start+w1-1, strand),
		      "...",
		      get_view(x, end-w2+1, end, strand), sep="")
	}
}


###
## Show header of output table
## NOTE: Based on function from Biostrings package
##
show_vframe_header <- function(iW, cols)
{
	cat(format("", width=iW)) # Print padding
	
	for (col in cols)
	{# Print column names
		cat(" ")
		cat(format(col$nm, width=col$width, justify="right"))
	}
	cat("\n")
}


###
## Show row of output table
## NOTE: Based on function from Biostrings package
##
show_vframe_line <- function(x, i, iW, cols)
{
	# Print triplex index
	cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"))
	
	colW <- 0 # Summary of all column width
	
	for (col in cols)
	{# Print column values
		cat(" ")
		value <- do.call(col$fn, list(x))[i]
		cat(do.call("format", c(value, col[3:length(col)], justify="right")))
		
		colW = colW + col$width
	}
	snippetW <- getOption("width") - iW - colW - length(cols) - 3
	cat(" [",
	    get_snippet(subject(x), start(x)[i], end(x)[i], snippetW, strand(x)[i]),
	    "]\n", sep="")
}


###
## Shot dots for not showed table rows
##
show_vframe_line_dots <- function(x, iW, cols)
{
	cat(format("...", width=iW, justify="right"))
	
	for (col in cols) {
		cat(" ")
		cat(format("...", width=col$width, justify="right"))
	}
	cat("\n")
}


###
## Show all output table rows
## NOTE: Based on function from Biostrings package
## 'half_nrow' must be >= 1
##
show_vframe <- function(x, half_nrow=9L)
{
	###
	## Column definitions
	## nm = Column header
	## fn = Function name to get row values
	## ... other parameters passed to format function
	##
	cols <- list(
		list(nm="start",  fn="start"),
		list(nm="width",  fn="width"),
		list(nm="score",  fn="score"),
		list(nm="pvalue", fn="pvalue", scientific=TRUE, digits=2),
		list(nm="ins",    fn="ins"),
		list(nm="type",   fn="type"),
		list(nm="s",      fn="strand")
	)
	
	i <- 1
	for (col in cols)
	{# Calculate column widths
		col_max <- max(do.call(col$fn, list(x)))
		col_maxstr <- do.call("format", c(col_max, col[3:length(col)]))
		cols[[i]] <- c(col, width=max(nchar(col_maxstr), nchar(col$nm)))
		i = i + 1
	}
	
	lx <- length(x)
	iW <- nchar(format(lx)) + 2 # Two extra for square brackets
	
	show_vframe_header(iW, cols)
	
	if (lx <= 2*half_nrow + 1)
	{# Show all
		for (i in seq_len(lx))
			show_vframe_line(x, i, iW, cols)
	}
	else
	{# Show first and last views
		for (i in 1:half_nrow)
			show_vframe_line(x, i, iW, cols)
		
		show_vframe_line_dots(x, iW, cols)
		
		for (i in (lx-half_nrow+1L):lx)
			show_vframe_line(x, i, iW, cols)
	}
}


###
## Show method, automatically invoked by print command
## NOTE: Based of function from Biostrings package
##
setMethod("show", "TriplexViews", function(object)
{
	subject <- subject(object)
	lsub <- length(subject)
	
	cat("  Triplex views on a ", lsub, "-letter ", class(subject),
			" subject", sep="")
	cat("\nsubject:", get_snippet(subject, 1, lsub, getOption("width") - 9, "+"))
	cat("\ntriplexes:")
	
	if (length(object) == 0) {
		cat(" NONE\n")
	}
	else {
		cat("\n")
		show_vframe(object)
	}
})


###
## Coerce TriplexViews to DNAStringSet
##
setAs("TriplexViews", "DNAStringSet", function(from)
{
	s <- DNAStringSet(subject(from), start(from), end(from))
	rev_cmpl <- strand(from) == "-"
	s[rev_cmpl] <- reverseComplement(s[rev_cmpl])
	
	for (i in 1:length(s))
	{# Set proper names
		names(s)[i] <- paste(
			"start=", start(from)[i], ";",
			"end=", end(from)[i], ";",
			"score=", score(from)[i], ";",
			"pvalue=", format(pvalue(from)[i], scientific=TRUE, digits=2), ";",
			"ins=", ins(from)[i], ";",
			"type=", type(from)[i], ";",
			"strand=", strand(from)[i], sep=""
		)
	}
	return(s)
})


###
## Coerce TriplexViews to GRanges
##
setAs("TriplexViews", "GRanges", function(from)
{
	seqlen <- length(subject(from))
	names(seqlen) <- "chr1"
	
	GRanges(
		"chr1",
		IRanges(start(from), end(from)),
		strand(from),
		score = score(from),
		tritype = type(from),
		pvalue = format(pvalue(from), scientific=TRUE, digits=2),
		lstart = lstart(from),
		lend = lend(from),
		indels = ins(from),
		seqlengths = seqlen
	)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Always behaves like an endomorphism (i.e. ignores the 'drop' argument and
### behaves like if it was actually set to FALSE).
setMethod("[", "TriplexViews", function(x, i, j, ..., drop=TRUE)
{
	if (!missing(j) || length(list(...)) > 0L)
		stop("invalid subsetting")
	if (missing(i))
		i <- seq_len(length(x))
	
	tx_views <- TriplexViews(
		subject(x),
		start  = start(x)[i],
		end    = end(x)[i],
		score  = score(x)[i],
		pvalue = pvalue(x)[i],
		insdel = ins(x)[i],
		type   = type(x)[i],
		lstart = lstart(x)[i],
		lend   = lend(x)[i],
		strand = strand(x)[i],
		params = x@params,
		score_table = x@score_table,
		group_table = x@group_table
	)
	return(tx_views)
})


###
## Coerce to character vector
##
setMethod("as.character", "TriplexViews", function(x)
{
	s <- as(x, "DNAStringSet")
	as.character(s)
})


###
## Convert to printable string
##
setMethod("toString", "TriplexViews", function(x, ...)
{
	toString(as.character(x), ...)
})
