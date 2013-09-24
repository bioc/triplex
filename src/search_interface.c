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
#include "XVector_interface.h"
#include "Biostrings_interface.h"

#include "search_interface.h"
#include "search.h"
#include "libtriplex.h"
#include "dl_list.h"


/* Global Variable  */
t_dl_list dl_list, dl_list_arr[8];
int act_dl_list;


/**
 * Decode DNAString object
 * @see IRanges_interface, Biostrings_interface
 * @param dnaobject DNAString object
 * @param seq_type Sequence type
 * @return Decoded sequence structure
 */
seq_t decode_DNAString(SEXP dnaobject, int seq_type)
{
	// Extract char sequence from R object
	cachedCharSeq x = cache_XRaw(dnaobject);
	// Initialize structure for decoded string
	seq_t dna;
	dna.len = x.length;
	dna.seq = malloc((x.length + 1) * sizeof(char));
	if (dna.seq == NULL)
		error("Failed to allocate memory for decoded DNA string.");
	
	double stat[NBASES];
	for (int j = 0; j < NBASES; j++)
	// Initialize statistics
		stat[j] = 0;
	
	int i; char ch;
	
	for (i = 0; i < dna.len; i++)
	{
		ch = CHAR2NUKL[tolower(DNAdecode(x.seq[i]))];
		if (ch == INVALID_CHAR)
		{
			free(dna.seq);
			error("Unsupported symbol '%c' in input sequence.", dna.seq[i]);
		}
		else
		{
			dna.seq[i] = ch;
			if (ch < NBASES) stat[ch]++;
		}
	}
	
	dna.seq[i] = '\0';
	dna.type = ST_PR;
	
	for (int j = 0; j < NBASES; j++)
	// Normalize statistics
		stat[j] = stat[j] / dna.len;
	
	Rprintf("seq_type = %d\n", seq_type);
	Rprintf("stat = %g %g %g %g\n", stat[A], stat[C], stat[G], stat[T]);
	
	return dna;
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
 * Search triplexes in DNA sequence
 * NOTE .Call entry point
 * @param dnaobject DNAString object
 * @param type      Triplex type vector
 * @param rparams   Custom algorithm options
 * @param pbw       Progress bar width
 * @return List
 */
SEXP triplex_search(SEXP dnaobject, SEXP type, SEXP seq_type, SEXP rparams, SEXP pbw)
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
	
	// Set Lambda, Mu and Rn tables
	LAMBDA[ST_PR][0] = LAMBDA[ST_PR][1] = LAMBDA[ST_PR][2] = LAMBDA[ST_PR][3] = p[P_LAMBDA_1_P];
	LAMBDA[ST_PR][4] = LAMBDA[ST_PR][5] = LAMBDA[ST_PR][6] = LAMBDA[ST_PR][7] = p[P_LAMBDA_2_P];
	LAMBDA[ST_EU][0] = LAMBDA[ST_EU][1] = LAMBDA[ST_EU][2] = LAMBDA[ST_EU][3] = p[P_LAMBDA_1_E];
	LAMBDA[ST_EU][4] = LAMBDA[ST_EU][5] = LAMBDA[ST_EU][6] = LAMBDA[ST_EU][7] = p[P_LAMBDA_2_E];
	
	MI[ST_PR][0] = MI[ST_PR][1] = MI[ST_PR][2] = MI[ST_PR][3] = p[P_MI_1_P];
	MI[ST_PR][4] = MI[ST_PR][5] = MI[ST_PR][6] = MI[ST_PR][7] = p[P_MI_2_P];
	MI[ST_EU][0] = MI[ST_EU][1] = MI[ST_EU][2] = MI[ST_EU][3] = p[P_MI_1_E];
	MI[ST_EU][4] = MI[ST_EU][5] = MI[ST_EU][6] = MI[ST_EU][7] = p[P_MI_2_E];
	
	RN[ST_PR][0] = RN[ST_PR][1] = RN[ST_PR][2] = RN[ST_PR][3] = p[P_RN_1_P];
	RN[ST_PR][4] = RN[ST_PR][5] = RN[ST_PR][6] = RN[ST_PR][7] = p[P_RN_2_P];
	RN[ST_EU][0] = RN[ST_EU][1] = RN[ST_EU][2] = RN[ST_EU][3] = p[P_RN_1_E];
	RN[ST_EU][4] = RN[ST_EU][5] = RN[ST_EU][6] = RN[ST_EU][7] = p[P_RN_2_E];
	
	int *st = INTEGER(seq_type);
	
	int *t = INTEGER(type);
	seq_t dna = decode_DNAString(dnaobject, st[0]);
	intv_t *chunk = get_chunks(dna);
	
	for (int i = 0; i < NUM_TRI_TYPES; i++)
		dl_list_init(&dl_list_arr[i], p[P_MAX_LEN]+p[P_MAX_LOOP]); // FIXME
	
	for (int i = 0; i < LENGTH(type); i++)
	{// Call original main function for all specified vector types
		act_dl_list = i;
		
		params.tri_type = t[i];
		main_search(dna, chunk, &params, &pen, *INTEGER(pbw));
		dl_list_group_filter(&dl_list_arr[i]);
	}
	
	dl_list_merge_sort(dl_list_arr, &dl_list, NUM_TRI_TYPES);
	list = export_results(&dl_list);
	dl_list_free(&dl_list);
	
	free(dna.seq);
	free_intv(chunk);
	
	return list;
}
