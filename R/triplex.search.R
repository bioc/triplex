###
## Triplex search R interface
##
## Author: Jiri Hon
## Date: 2012/10/15
## Package: triplex
##

###
## Check if type vector is valid
##
validate_type <- function(type)
{
	if (mode(type) != "numeric")
		stop("Triplex type must be a numeric vector, see ?triplex.search")
	
	type <- unique(as.integer(type))
		
	if (length(type) > 8 ||
	    min(type) < 0 ||
	    max(type) > 7)
	{
		stop("Invalid triplex type vector, see ?triplex.search")
	}
	return(type)
}

###
## Convert storage mode to double if possible
##
to_double <- function(val)
{
	if (mode(val) != "numeric")
		stop("Invalid parameter, use numeric values please.")
	
	return(as.double(val))
}

###
## Search input sequence for triplexes
##
triplex.search <- function(
	dna,
	type       = 0:7,
	min_score  = 10,
	p_value    = 1,
	min_len    = 6,
	max_len    = 25,
	min_loop   = 3,
	max_loop   = 10,
	lambda_1   = 0.67,
	lambda_2   = 0.71,
	mi_1       = 6.05,
	mi_2       = 5.88,
	dtwist_pen = 7,
	ins_pen    = 9,
	iso_pen    = 5,
	iso_bonus  = 0,
	mis_pen    = 7)
{
	if (class(dna) != "DNAString")
		stop("Input sequence must be DNAString object.")
	
	if (min_loop < 1)
		stop("Can not search triplexes whithout a loop.")
	
	if (max_loop < min_loop)
		stop("max_loop option can not be lower than min_loop.")
	
	if (min_len < 1)
		stop("Can not search triplexes whithout a body.")
	
	if (max_len < min_len)
		stop("max_len option can not be lower than min_len.")
	
	# Check and convert parameters to double
	p <- double()
	p[MIN_SCORE]  = to_double(min_score)
	p[P_VALUE]    = to_double(p_value)
	p[MIN_LEN]    = to_double(min_len)
	p[MAX_LEN]    = to_double(max_len)
	p[MIN_LOOP]   = to_double(min_loop)
	p[MAX_LOOP]   = to_double(max_loop)
	p[LAMBDA_1]   = to_double(lambda_1)
	p[LAMBDA_2]   = to_double(lambda_2)
	p[MI_1]       = to_double(mi_1)
	p[MI_2]       = to_double(mi_2)
	p[DTWIST_PEN] = to_double(dtwist_pen)
	p[INS_PEN]    = to_double(ins_pen)
	p[ISO_PEN]    = to_double(iso_pen)
	p[ISO_BONUS]  = to_double(iso_bonus)
	p[MIS_PEN]    = to_double(mis_pen)
	
	type <- validate_type(type)
	
	txs <- .Call("triplex_search", dna, type, p, as.integer(getOption("width")))
	
	strand <- txs[[T_STRAND]]
	s <- character()
	s[strand == 0] <- "+"
	s[strand == 1] <- "-"
	
	tx_views <- TriplexViews(
		dna,
		start  = txs[[T_START]],
		end    = txs[[T_END]],
		score  = txs[[T_SCORE]],
		pvalue = txs[[T_P_VALUE]],
		insdel = txs[[T_INSDEL]],
		type   = txs[[T_TYPE]],
		lstart = txs[[T_L_START]],
		lend   = txs[[T_L_END]],
		strand = s,
		params = p
	)
	return(tx_views)
}
