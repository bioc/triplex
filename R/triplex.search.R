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
## Validate sequence type
##
validate_seq_type <- function(seq_type)
{
	# Validate sequence type option
	if (seq_type == 'prokaryotic')
		seq_type = as.integer(ST_PROKARYOTIC)
	else if (seq_type == 'eukaryotic')
		seq_type = as.integer(ST_EUKARYOTIC)
	else
		error("Invalid sequence type. Valid options are prokaryotic or eukaryotic.")
}

###
## Validate score or group table
##
validate_table <- function(table, type)
{
	if (table == 'default')
	{# Return default tables
		if (type == 'score')
			return(triplex.score.table())
		else if (type == 'group')
			return(triplex.group.table())
	}
	
	if (is.null(table$par) || is.null(table$apar))
		error("Incomplete custom table.")
	
	if (!setequal(dim(table$par), dim(table$apar)) || !setequal(dim(table$par), c(NBASES, NBASES)))
		error("Invalid size of custom table.")
	
	table$par <- as.integer(table$par)
	table$apar <- as.integer(table$apar)
	
	return(table)
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
		return(value)
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
	type        = 0:7,
	min_score   = 10,
	p_value     = 1,
	min_len     = 6,
	max_len     = 25,
	min_loop    = 3,
	max_loop    = 10,
	seq_type    = 'eukaryotic',
	score_table = 'default',
	group_table = 'default',
	lambda_par  = 'default',
	lambda_apar = 'default',
	mu_par      = 'default',
	mu_apar     = 'default',
	rn_par      = 'default',
	rn_apar     = 'default',
	dtwist_pen  = 'default', #7,
	ins_pen     = 'default', #9,
	iso_pen     = 'default', #5,
	iso_bonus   = 'default', #0,
	mis_pen     = 'default') #7)
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
			'Default penalization options were changed. You should consider recalculation
			of P-value constants (mi, lambda, rn) . For details see vignette.')
	
	# Set default penalizations
	dtwist_pen = pen_set(dtwist_pen, 7)
	ins_pen    = pen_set(ins_pen, 9)
	iso_pen    = pen_set(iso_pen, 5)
	iso_bonus  = pen_set(iso_bonus, 0)
	mis_pen    = pen_set(mis_pen, 7)
	
	###
	## Set default P-value constants
	##
	## Triplexfit output as reported by Matej Lexa on 2013-04-30:
	##
	## Sequence length: 4800000
	##
	## Ecoli_par:  mu = 7.4805, lambda = 0.8892, hits = 194641, rn = 0.0406
	## Ecoli_apar: mu = 7.6569, lambda = 0.8092, hits = 130954, rn = 0.0273
	## Human_par:  mu = 7.5835, lambda = 0.8433, hits = 145915, rn = 0.0304
	## Human_apar: mu = 7.9611, lambda = 0.6910, hits = 194391, rn = 0.0405
	
	lambda_par = pconst_set(lambda_par, 0.8892, 0.8433)
	lambda_apar = pconst_set(lambda_apar, 0.8092, 0.6910)
	mu_par = pconst_set(mu_par, 7.4805, 7.5835)
	mu_apar = pconst_set(mu_apar, 7.6569, 7.9611)
	rn_par = pconst_set(rn_par, 0.0406, 0.0304)
	rn_apar = pconst_set(rn_apar, 0.0273, 0.0405)
	
	# Check and convert parameters to double
	p <- double()
	p[MIN_SCORE]     = to_double(min_score)
	p[P_VALUE]       = to_double(p_value)
	p[MIN_LEN]       = to_double(min_len)
	p[MAX_LEN]       = to_double(max_len)
	p[MIN_LOOP]      = to_double(min_loop)
	p[MAX_LOOP]      = to_double(max_loop)
	p[LAMBDA_PAR_P]  = to_double(lambda_par[1])
	p[LAMBDA_PAR_E]  = to_double(lambda_par[2])
	p[LAMBDA_APAR_P] = to_double(lambda_apar[1])
	p[LAMBDA_APAR_E] = to_double(lambda_apar[2])
	p[MU_PAR_P]      = to_double(mu_par[1])
	p[MU_PAR_E]      = to_double(mu_par[2])
	p[MU_APAR_P]     = to_double(mu_apar[1])
	p[MU_APAR_E]     = to_double(mu_apar[2])
	p[RN_PAR_P]      = to_double(rn_par[1])
	p[RN_PAR_E]      = to_double(rn_par[2])
	p[RN_APAR_P]     = to_double(rn_apar[1])
	p[RN_APAR_E]     = to_double(rn_apar[2])
	p[DTWIST_PEN]    = to_double(dtwist_pen)
	p[INS_PEN]       = to_double(ins_pen)
	p[ISO_PEN]       = to_double(iso_pen)
	p[ISO_BONUS]     = to_double(iso_bonus)
	p[MIS_PEN]       = to_double(mis_pen)
	
	type <- validate_type(type)
	seq_type <- validate_seq_type(seq_type)
	score_table <- validate_table(score_table, 'score')
	group_table <- validate_table(group_table, 'group')
	
	txs <- .Call(
		"triplex_search", dna, type, seq_type, p,
		score_table$par, score_table$apar,
		group_table$par, group_table$apar,
		as.integer(getOption("width")))
	
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
