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
## Set penalization value
##
pen_set <- function(value, default)
{
	if (value == 'default')
		return(default)
	else
		retrun(value)
}

###
## Set P-value constant
##
pconst_set <- function(value, prokaryotic, eukaryotic)
{
	if (value == 'default')
		return(c(prokaryotic, eukaryotic))
	else
		return(c(value, value))
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
	seq_type   = 'auto',
	lambda_1   = 'default',
	lambda_2   = 'default',
	mi_1       = 'default',
	mi_2       = 'default',
	rn_1       = 'default',
	rn_2       = 'default',
	dtwist_pen = 'default', #7,
	ins_pen    = 'default', #9,
	iso_pen    = 'default', #5,
	iso_bonus  = 'default', #0,
	mis_pen    = 'default', #7)
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
	
	if (dtwist_pen != 'default' || ins_pen != 'default' ||
		 iso_pen != 'default' || iso_bonus != 'default' ||
		 mis_pen != 'default')
		warning(
			'Default penalization options were changed, so you should consider P-value
			constants (mi, lambda, rn) recalculation. For details see vignette.')
	
	# Set default penalizations
	dtwist_pen = pen_set(dtwist_pen, 7)
	ins_pen    = pen_set(ins_pen, 9)
	iso_pen    = pen_set(iso_pen, 5)
	iso_bonus  = pen_set(iso_bonus, 0)
	mis_pen    = pen_set(mis_pen, 7)
	
	# Validate sequence type option
	if (seq_type == 'auto')
		seq_type = S_AUTO
	else if (seq_type == 'prokaryotic')
		seq_type = S_PROKARYOTIC
	else if (seq_type == 'eukaryotic')
		seq_type = S_EUKARYOTIC
	else
		error("Invalid sequence type. Valid options: auto, prokaryotic, eukaryotic")
	
	# Set default P-value constants
	#
	# Triplexfit output as reported by Matej Lexa on 2013-04-30:
	#
	# Sequence length: 4800000
	#
	# Ecoli_par:  mu = 7.4805, lambda = 0.8892, hits = 194641, rn = 0.0406
	# Ecoli_apar: mu = 7.6569, lambda = 0.8092, hits = 130954, rn = 0.0273
	# Human_par:  mu = 7.5835, lambda = 0.8433, hits = 145915, rn = 0.0304
	# Human_apar: mu = 7.9611, lambda = 0.6910, hits = 194391, rn = 0.0405
	
	lambda_1 = pconst_set(lambda_1, 0.8892, 0.8433)
	lambda_2 = pconst_set(lambda_2, 0.8092, 0.6910)
	mi_1 = pconst_set(mi_1, 7.4805, 7.5835)
	mi_2 = pconst_set(mi_2, 7.6569, 7.9611)
	rn_1 = pconst_set(rn_1, 0.0406, 0.0304)
	rn_2 = pconst_set(rn_2, 0.0273, 0.0405)
	
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
