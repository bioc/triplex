/**
 * Triplex package
 * C interface to search algorithm
 *
 * @author  Jiri Hon, Tomas Martinek
 * @date    2012/10/15
 * @file    search_interface.c
 * @package triplex
 */

#include <ctype.h>

#include "IRanges_interface.h"
#include "Biostrings_interface.h"

#include "search_interface.h"
#include "search.h"
#include "libtriplex.h"
#include "dl_list.h"

#define ASCII_LOW 128
#define INVALID_CHAR -1

/* Global Variable  */
t_dl_list dl_list, dl_list_arr[8];
int act_dl_list;

/* Translation table from ascii symbols to
 * internal representation of DNA bases */
int CHAR2NUKL[ASCII_LOW];


/**
 * Decode DNAString object
 * @see IRanges_interface, Biostrings_interface
 * @param dnaobject DNAString object
 * @return Decoded sequence structure
 */
seq_t decode_DNAString(SEXP dnaobject)
{
	// Extract char sequence from R object
	cachedCharSeq x = cache_XRaw(dnaobject);
	// Initialize structure for decoded string
	seq_t dna;
	dna.seq = Calloc(x.length + 1, char);
	dna.len = x.length;
	
	int i;
	for (i = 0; i < dna.len; i++)
		dna.seq[i] = tolower(DNAdecode(x.seq[i]));
	
	dna.seq[i] = '\0';
	return dna;
}


/**
 * Encode bases to enum values A,C,G,T
 * @param seq DNA sequence
 */
void encode_bases(char *seq)
{
	int ch;
	for (int i = 0; seq[i] != '\0'; i++)
	{
		ch = CHAR2NUKL[(unsigned) seq[i]];
		if (ch == INVALID_CHAR)
			error("Unsupported symbol '%c' in input sequence.", seq[i]);
		else
			seq[i] = (char) ch;
	}
}


/**
 * Create list element
 * @param list List object
 * @param idx  Element index
 * @param type Vector type
 * @param len  Vector length
 * @return Vector object
 */
SEXP create_list_elt(SEXP list, int idx, SEXPTYPE type, int len)
{
	SEXP vector;
	PROTECT(vector = allocVector(type, len));
	SET_VECTOR_ELT(list, idx, vector);
	UNPROTECT(1);
	return vector;
}


/**
 * Export results in list object
 * @param dl_list List pointer
 * @return List object
 */
SEXP export_results(t_dl_list *dl_list)
{
	t_dl_node *pointer;
	SEXP list;
	PROTECT(list = allocVector(VECSXP, 9));
	
	int *start = INTEGER(create_list_elt(list, 0, INTSXP, dl_list->size));
	int *end = INTEGER(create_list_elt(list, 1, INTSXP, dl_list->size));
	int *score = INTEGER(create_list_elt(list, 2, INTSXP, dl_list->size));
	double *pvalue = REAL(create_list_elt(list, 3, REALSXP, dl_list->size));
	int *insdel = INTEGER(create_list_elt(list, 4, INTSXP, dl_list->size));
	int *type = INTEGER(create_list_elt(list, 5, INTSXP, dl_list->size));
	int *lstart = INTEGER(create_list_elt(list, 6, INTSXP, dl_list->size));
	int *lend = INTEGER(create_list_elt(list, 7, INTSXP, dl_list->size));
	int *strand = INTEGER(create_list_elt(list, 8, INTSXP, dl_list->size));
	
	pointer = (dl_list->first)->next;
	
	/* printf("List size: %d\n", dl_list->size); */
	
	for (int i = 0; i < dl_list->size; i++)
	{
		/* printf("(%d,%d)", pointer->data.start, pointer->data.end); */
		start[i] = pointer->data.start;
		end[i] = pointer->data.end;
		score[i] = pointer->data.score;
		pvalue[i] = pointer->data.pvalue;
		insdel[i] = pointer->data.insdel;
		type[i] = pointer->data.type;
		lstart[i] = pointer->data.lstart;
		lend[i] = pointer->data.lend;
		strand[i] = pointer->data.strand;

                /* Rprintf("Triplex: %d: %d, %d, %d, %d\n", score[i], start[i], lstart[i], lend[i], end[i]); */
		
		pointer = pointer->next;
	}
	
	/* printf("\nEnd\n"); */
	
	UNPROTECT(1);
	return list;
}


/**
 * Save result
 * @param start  Start of triplex
 * @param end    End of triplex
 * @param score  Triplex score
 * @param pvalue P-value
 * @param insdel Number of insertions/deletions
 * @param type   Triplex type
 * @param lstart Loop start
 * @param lend   Loop end
 */
void save_result(
	int start, int end, int score, double pvalue, int insdel,
	int type,  int lstart, int lend, int strand)
{
	t_dl_data data =
	{
		.start = start,
		.end = end,
		.score = score,
		.pvalue = pvalue,
		.insdel = insdel,
		.type = type,
		.lstart = lstart,
		.lend = lend,
		.strand = strand 
	};
	dl_list_insert(&dl_list_arr[act_dl_list], data);
}


/**
 * Initialize global CHAR2NUKL translation table
 * Called from @see R_init_triplex
 */
void init_CHAR2NUKL_table()
{
	for (int i = 0; i < ASCII_LOW; i++)
		CHAR2NUKL[i] = INVALID_CHAR;
	
	CHAR2NUKL['r'] = 'r';
	CHAR2NUKL['m'] = 'm';
	CHAR2NUKL['w'] = 'w';
	CHAR2NUKL['d'] = 'd';
	CHAR2NUKL['v'] = 'v';
	CHAR2NUKL['h'] = 'h';
	CHAR2NUKL['b'] = 'b';
	CHAR2NUKL['s'] = 's';
	CHAR2NUKL['y'] = 'y';
	CHAR2NUKL['k'] = 'k';
	CHAR2NUKL['a'] = A;
	CHAR2NUKL['c'] = C;
	CHAR2NUKL['g'] = G;
	CHAR2NUKL['t'] = T;
}


/**
 * Search triplexes in DNA sequence
 * .Call entry point
 * @param dnaobject DNAString object
 * @param type      Triplex type vector
 * @param rparams   Custom algorithm options
 * @param pbw       Progress bar width
 * @return List
 */
SEXP triplex_search(SEXP dnaobject, SEXP type, SEXP rparams, SEXP pbw)
{
   SEXP list;
	
	double *p = REAL(rparams);
	
	t_params params =
	{// Set params
		.tri_type = -1,
		.min_score = p[P_MIN_SCORE],
		.p_val = p[P_P_VALUE],
		.min_len = p[P_MIN_LEN],
		.max_len = p[P_MAX_LEN],
		.min_loop = p[P_MIN_LOOP],
		.max_loop = p[P_MAX_LOOP],
		.max_chunk_size = MAX_CHUNK_SIZE
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
	
	int *t = INTEGER(type);
	seq_t dna;
	
	for (int i = 0; i<8; i++)
		dl_list_init(&dl_list_arr[i], p[P_MAX_LEN]+p[P_MAX_LOOP]); // FIXME
	
	for (int i = 0; i < LENGTH(type); i++)
	{// Call original main function for all specified vector types
		act_dl_list = i;

		dna = decode_DNAString(dnaobject);
		
		params.tri_type = t[i];
		main_search(dna, params, pen, *INTEGER(pbw));
		dl_list_group_filter(&dl_list_arr[i]);
		
		Free(dna.seq);
	}
	dl_list_merge_sort(dl_list_arr, &dl_list, 8);
	list = export_results(&dl_list);
	dl_list_free(&dl_list);
	
	return list;
}
