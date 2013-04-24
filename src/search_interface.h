/**
 * Triplex package
 * Header file for C interface to search algorithm
 *
 * @author  Jiri Hon
 * @date    2012/10/15
 * @file    search_interface.h
 * @package triplex
 */

#ifndef SEARCH_INTERFACE_H
#define SEARCH_INTERFACE_H

#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

#include "libtriplex.h"


typedef enum
{// Enumeration for decoding params from R
	P_MIN_SCORE = 0,
	P_P_VALUE,
	P_MIN_LEN,
	P_MAX_LEN,
	P_MIN_LOOP,
	P_MAX_LOOP,
	P_LAMBDA_1,
	P_LAMBDA_2,
	P_MI_1,
	P_MI_2,
	P_DTWIST_PEN,
	P_INS_PEN,
	P_ISO_PEN,
	P_ISO_BONUS,
	P_MIS_PEN
} rparams_t;


SEXP triplex_search(SEXP dnaobject, SEXP type, SEXP params, SEXP pbw);
seq_t decode_DNAString(SEXP dnaobject);
void save_result(
	int start, int end,    int score, double pvalue, int insdel,
	int type,  int lstart, int lend, int strand
);

#endif // SEARCH_INTERFACE_H
