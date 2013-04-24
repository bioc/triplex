/**
 * Triplex package
 * C interface for calling alignment function from R
 *
 * @author  Jiri Hon
 * @date    2012/10/15
 * @file    align_interface.c
 * @package triplex
 */

#include "align.h"
#include "align_interface.h"
#include "search.h"
#include "search_interface.h"

seq_t tx_align;
int tx_align_pos;


/**
 * Print triplex alignment into global string
 * @param fmt Format string
 * @param ... Other parameters
 */
void Aprintf(char ch)
{
	if (tx_align_pos < tx_align.len)
		tx_align.seq[tx_align_pos++] = ch;
}


/**
 * Align triplex sequence
 * .Call entry point
 * @param seq     Character vector with sequence
 * @param type    Triplex type
 * @param rparams Params of search algorithm
 * @return aligned sequence for visualization
 */
SEXP triplex_align(SEXP seq, SEXP type, SEXP rparams)
{
	double *p = REAL(rparams);
	double *t = REAL(type);
	
	t_params params =
	{// Set params
		.tri_type = t[0],
		.min_score = p[P_MIN_SCORE],
		.p_val = p[P_P_VALUE],
		.min_len = p[P_MIN_LEN],
		.max_len = p[P_MAX_LEN],
		.min_loop = p[P_MIN_LOOP],
		.max_loop = p[P_MAX_LOOP]
	};
	t_penalization pen =
	{// Set penalizations 
		.dtwist = p[P_DTWIST_PEN],
		.insertion = p[P_INS_PEN],
		.iso_change = p[P_ISO_PEN],
		.iso_stay = p[P_ISO_BONUS],
		.mismatch = p[P_MIS_PEN] 
	};
	// Set Lambda and Mi tables
	LAMBDA[0] = LAMBDA[1] = LAMBDA[4] = LAMBDA[5] = p[P_LAMBDA_1];
	LAMBDA[2] = LAMBDA[3] = LAMBDA[6] = LAMBDA[7] = p[P_LAMBDA_2];
	MI[0] = MI[1] = MI[4] = MI[5] = p[P_MI_1];
	MI[2] = MI[3] = MI[6] = MI[7] = p[P_MI_2];
	
	seq_t dna = decode_DNAString(seq);
	
	// Ininitalize global alignment string
	tx_align.len = 2 * dna.len;
	tx_align.seq = Calloc(tx_align.len, char);
	tx_align_pos = 0;
	
	main_align(dna, params, pen);
	
	SEXP res;
	PROTECT(res = allocVector(STRSXP, 1));
	SET_STRING_ELT(res, 0, mkChar(tx_align.seq));
	UNPROTECT(1);
	
	Free(dna.seq);
	Free(tx_align.seq);
	
	return res;
}
