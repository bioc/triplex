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
 * Fill scoring and group tables based on user input from R
 * @param score_par Score table for parallels
 * @param score_apar Score table for antiparallels
 * @param group_par Isogroup table for parallels
 * @param group_apar Isogroup table for antiparallels
 */
void score_group_tables_fill(
	int score_par[NBASES][NBASES], int score_apar[NBASES][NBASES],
	int group_par[NBASES][NBASES], int group_apar[NBASES][NBASES])
{
	int i, j;
	
	/* Type 0 */
	for (i = 0; i < NBASES; i++)
		for (j = 0; j < NBASES; j++)
		{
			TAB_SCORE[0][i][j] = score_par[i][COMP[j]];
			TAB_GROUP[0][i][j] = group_par[i][COMP[j]];
		}
	
	/* Type 1 */
	for (i = 0; i < NBASES; i++)
		for (j = 0; j < NBASES; j++)
		{
			TAB_SCORE[1][i][j] = score_par[j][COMP[i]];
			TAB_GROUP[1][i][j] = group_par[j][COMP[i]];
		}
	
	/* Type 2 */
	for (i = 0; i < NBASES; i++)
		for (j = 0; j < NBASES; j++)
		{
			TAB_SCORE[2][i][j] = score_par[COMP[i]][j];
			TAB_GROUP[2][i][j] = group_par[COMP[i]][j];
		}
	
	/* Type 3 */
	for (i = 0; i < NBASES; i++)
		for (j = 0; j < NBASES; j++)
		{
			TAB_SCORE[3][i][j] = score_par[COMP[j]][i];
			TAB_GROUP[3][i][j] = group_par[COMP[j]][i];
		}
	
	/* Type 4 */
	for (i = 0; i < NBASES; i++)
		for (j = 0; j < NBASES; j++)
		{
			TAB_SCORE[4][i][j] = score_apar[COMP[i]][COMP[j]];
			TAB_GROUP[4][i][j] = group_apar[COMP[i]][COMP[j]];
		}
	
	/* Type 5 */
	for (i = 0; i < NBASES; i++)
		for (j = 0; j < NBASES; j++)
		{
			TAB_SCORE[5][i][j] = score_apar[COMP[j]][COMP[i]];
			TAB_GROUP[5][i][j] = group_apar[COMP[j]][COMP[i]];
		}
	
	/* Type 6 */
	for (i = 0; i < NBASES; i++)
		for (j = 0; j < NBASES; j++)
		{
			TAB_SCORE[6][i][j] = score_apar[i][j];
			TAB_GROUP[6][i][j] = group_apar[i][j];
		}
	
	/* Type 7 */
	for (i = 0; i < NBASES; i++)
		for (j = 0; j < NBASES; j++)
		{
			TAB_SCORE[7][i][j] = score_apar[j][i];
			TAB_GROUP[7][i][j] = group_apar[j][i];
		}
}


/**
 * Save integer R-matrix saved by columns into integer C-array saved by rows
 * @param r_matrix R-matrix data saved by columns
 * @param c_matrix Output C-array saved by rows
 * @param nrow Number of rows
 * @param ncol Number of columns
 */
void col_matrix_to_row_array(int *r_matrix, int *c_array, int nrow, int ncol)
{
	for (int c = 0; c < ncol; c++)
		for (int r = 0; r < nrow; r++)
			c_array[r*ncol + c] = r_matrix[c*nrow + r];
}


/**
 * Print integer table
 * @param table Table data pointer
 * @param nrow Number of rows
 * @param ncol Number of columns
 */
void print_table(int *table, int nrow, int ncol)
{
	for (int i = 0; i < nrow; i++)
	{
		for (int j = 0; j < ncol; j++)
			Rprintf("%3d ", table[i*nrow + j]);
		
		Rprintf("\n");
	}
}


/**
 * Check two sets of score and group tables for identity
 * @param score_old Old score table
 * @param score_new New score table
 * @param group_old Old isogroup table
 * @param group_new New isogroup table
 */
void compare_score_group_tables(
	int score_old[NUM_TRI_TYPES][NBASES][NBASES],
	int score_new[NUM_TRI_TYPES][NBASES][NBASES],
	int group_old[NUM_TRI_TYPES][NBASES][NBASES],
	int group_new[NUM_TRI_TYPES][NBASES][NBASES])
{
	int score_marks[NUM_TRI_TYPES], group_marks[NUM_TRI_TYPES];
	memset(score_marks, 0, NUM_TRI_TYPES*sizeof(int));
	memset(group_marks, 0, NUM_TRI_TYPES*sizeof(int));
	
	for (int t = 0; t < NUM_TRI_TYPES; t++)
		for (int i = 0; i < NBASES; i++)
			for (int j = 0; j < NBASES; j++)
			{
				if (score_old[t][i][j] != score_new[t][i][j])
					score_marks[t] = 1;
				if (group_old[t][i][j] != group_new[t][i][j])
					group_marks[t] = 1;
			}
	
	int score_diff = 0, group_diff = 0;
	
	for (int t = 0; t < NUM_TRI_TYPES; t++)
	{
		if (score_marks[t])
		{
			score_diff = 1;
			Rprintf("Old score for %d:\n", t);
			print_table((int *) score_old[t], NBASES, NBASES);
			Rprintf("New score for %d:\n", t);
			print_table((int *) score_new[t], NBASES, NBASES);
		}
		if (group_marks[t])
		{
			group_diff = 1;
			Rprintf("Old group for %d:\n", t);
			print_table((int *) group_old[t], NBASES, NBASES);
			Rprintf("New group for %d:\n", t);
			print_table((int *) group_new[t], NBASES, NBASES);
		}
	}
	if (!score_diff && !group_diff)
		Rprintf("Success, no difference.\n");
}

#ifndef NDEBUG
/* Lookup tables */
int tab_score[NUM_TRI_TYPES][4][4] = {
 /*
  * Parallel
  */

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {-9, -9, -9, -9},
 /*   C   */    {-9,  2, -9, -9},
 /*   G   */    { 2,  1, -9, -9},
 /*   T   */    {-9,  1,  1,  2} },

 /* [W][H] */ /*  A   C   G   T */
 /*   A   */  { {-9, -9,  2, -9},
 /*   C   */    {-9,  2,  1,  1},
 /*   G   */    {-9, -9, -9,  1},
 /*   T   */    {-9, -9, -9,  2} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { { 2,  1,  1, -9},
 /*   C   */    {-9, -9,  1,  2},
 /*   G   */    {-9, -9,  2, -9},
 /*   T   */    {-9, -9, -9, -9} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { { 2, -9, -9, -9},
 /*   C   */    { 1, -9, -9, -9},
 /*   G   */    { 1,  1,  2, -9},
 /*   T   */    {-9,  2, -9, -9} },

 /*
  * Antiparallel
  */

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {-9, -9,  1,  2},
 /*   C   */    {-9,  2, -9, -9},
 /*   G   */    {-9, -9, -9,  1},
 /*   T   */    {-9,  1, -9,  2} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {-9, -9, -9, -9},
 /*   C   */    {-9,  2, -9,  1},
 /*   G   */    { 1, -9, -9, -9},
 /*   T   */    { 2, -9,  1,  2} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { { 2, -9,  1, -9},
 /*   C   */    { 1, -9, -9, -9},
 /*   G   */    {-9, -9,  2, -9},
 /*   T   */    { 2,  1, -9, -9} },

 /* [W][H] */ /*  A   C   G   T */
 /*   A   */  { { 2,  1, -9,  2},
 /*   C   */    {-9, -9, -9,  1},
 /*   G   */    { 1, -9,  2, -9},
 /*   T   */    {-9, -9, -9, -9} },

};

int tab_bound[NUM_TRI_TYPES][4][4] = {
 /*
  * Parallel
  */

 /* [H][W] */ /*  A   C   G   A */
 /*   A   */  { { 0,  0,  0,  0},
 /*   C   */    { 0,  1,  0,  0},
 /*   G   */    { 2,  2,  0,  0},
 /*   T   */    { 0,  1,  2,  1} },

 /* [W][H] */ /*  A   C   G   T */
 /*   A   */  { { 0,  0,  2,  0},
 /*   C   */    { 0,  1,  2,  1},
 /*   G   */    { 0,  0,  0,  2},
 /*   T   */    { 0,  0,  0,  1} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { { 1,  2,  1,  0},
 /*   C   */    { 0,  0,  2,  2},
 /*   G   */    { 0,  0,  1,  0},
 /*   T   */    { 0,  0,  0,  0} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { { 1,  0,  0,  0},
 /*   C   */    { 2,  0,  0,  0},
 /*   G   */    { 1,  2,  1,  0},
 /*   T   */    { 0,  2,  0,  0} },

 /*
  * Antiparallel
  */

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { { 0,  0,  3,  1},
 /*   C   */    { 0,  3,  0,  0},
 /*   G   */    { 0,  0,  0,  2},
 /*   T   */    { 0,  2,  0,  1} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { { 0,  0,  0,  0},
 /*   C   */    { 0,  3,  0,  2},
 /*   G   */    { 3,  0,  0,  0},
 /*   T   */    { 1,  0,  2,  1} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { { 1,  0,  2,  0},
 /*   C   */    { 2,  0,  0,  0},
 /*   G   */    { 0,  0,  3,  0},
 /*   T   */    { 1,  3,  0,  0} },

 /* [W][H] */ /*  A   C   G   T */
 /*   A   */  { { 1,  2,  0,  1},
 /*   C   */    { 0,  0,  0,  3},
 /*   G   */    { 2,  0,  3,  0},
 /*   T   */    { 0,  0,  0,  0} },

};
#endif

/**
 * Set score and group tables
 */
void set_score_group_tables(
	int *st_par, int *st_apar, int *gt_par, int *gt_apar)
{
	int score_par[NBASES][NBASES], score_apar[NBASES][NBASES],
	group_par[NBASES][NBASES], group_apar[NBASES][NBASES];
	
	/* Convert between R and C representation of matrix */
	col_matrix_to_row_array(st_par, (int *) score_par, NBASES, NBASES);
	col_matrix_to_row_array(st_apar, (int *) score_apar, NBASES, NBASES);
	col_matrix_to_row_array(gt_par, (int *) group_par, NBASES, NBASES);
	col_matrix_to_row_array(gt_apar, (int *) group_apar, NBASES, NBASES);
	
#ifndef NDEBUG
	int score_backup[NUM_TRI_TYPES][NBASES][NBASES];
	memcpy(score_backup, TAB_SCORE, NUM_TRI_TYPES*NBASES*NBASES*sizeof(int));
	int group_backup[NUM_TRI_TYPES][NBASES][NBASES];
	memcpy(group_backup, TAB_GROUP, NUM_TRI_TYPES*NBASES*NBASES*sizeof(int));
#endif
	
	score_group_tables_fill(score_par, score_apar, group_par, group_apar);

#ifndef NDEBUG
	compare_score_group_tables(tab_score, TAB_SCORE, tab_bound, TAB_GROUP);
	compare_score_group_tables(score_backup, TAB_SCORE, group_backup, TAB_GROUP);
#endif
}


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
	Chars_holder x = hold_XRaw(dnaobject);
	// Initialize structure for decoded string
	seq_t dna;
	dna.len = x.length;
	dna.seq = malloc((x.length + 1) * sizeof(char));
	if (dna.seq == NULL)
		error("Failed to allocate memory for decoded DNA string.");
	
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
		}
	}
	
	dna.seq[i] = '\0';
	dna.type = seq_type;
	
	return dna;
}


/**
 * Create list element
 * @param list List object
 * @param idx  Element index
 * @param type Vector type
 * @param len Vector length
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
 * Set lambda, mu and rn static tables
 * @param p Parameter vector
 */
void set_lambda_mu_rn_tables(double *p)
{
	LAMBDA[ST_PR][0] = LAMBDA[ST_PR][1] = LAMBDA[ST_PR][2] = LAMBDA[ST_PR][3] = p[P_LAMBDA_PAR_P];
	LAMBDA[ST_PR][4] = LAMBDA[ST_PR][5] = LAMBDA[ST_PR][6] = LAMBDA[ST_PR][7] = p[P_LAMBDA_APAR_P];
	LAMBDA[ST_EU][0] = LAMBDA[ST_EU][1] = LAMBDA[ST_EU][2] = LAMBDA[ST_EU][3] = p[P_LAMBDA_PAR_E];
	LAMBDA[ST_EU][4] = LAMBDA[ST_EU][5] = LAMBDA[ST_EU][6] = LAMBDA[ST_EU][7] = p[P_LAMBDA_APAR_E];
	
	MI[ST_PR][0] = MI[ST_PR][1] = MI[ST_PR][2] = MI[ST_PR][3] = p[P_MI_PAR_P];
	MI[ST_PR][4] = MI[ST_PR][5] = MI[ST_PR][6] = MI[ST_PR][7] = p[P_MI_APAR_P];
	MI[ST_EU][0] = MI[ST_EU][1] = MI[ST_EU][2] = MI[ST_EU][3] = p[P_MI_PAR_E];
	MI[ST_EU][4] = MI[ST_EU][5] = MI[ST_EU][6] = MI[ST_EU][7] = p[P_MI_APAR_E];
	
	RN[ST_PR][0] = RN[ST_PR][1] = RN[ST_PR][2] = RN[ST_PR][3] = p[P_RN_PAR_P];
	RN[ST_PR][4] = RN[ST_PR][5] = RN[ST_PR][6] = RN[ST_PR][7] = p[P_RN_APAR_P];
	RN[ST_EU][0] = RN[ST_EU][1] = RN[ST_EU][2] = RN[ST_EU][3] = p[P_RN_PAR_E];
	RN[ST_EU][4] = RN[ST_EU][5] = RN[ST_EU][6] = RN[ST_EU][7] = p[P_RN_APAR_E];
}


/**
 * Search triplexes in DNA sequence
 * NOTE .Call entry point
 * @param dnaobject DNAString object
 * @param type      Triplex type vector
 * @param rparams   Custom algorithm options
 * @param st_par Score table for parallel triplexes
 * @param st_apar Score table for antiparallel triplexes
 * @param gt_par Isogroup table for parallel triplexes
 * @param gt_apar Isogroup table for antiparallel triplexes
 * @param pbw       Progress bar width
 * @return List
 */
SEXP triplex_search(
	SEXP dnaobject, SEXP type, SEXP seq_type, SEXP rparams,
	SEXP st_par, SEXP st_apar, SEXP gt_par, SEXP gt_apar,
	SEXP pbw)
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
	
	int *st = INTEGER(seq_type);
	int *t = INTEGER(type);
	
	set_lambda_mu_rn_tables(p);
	set_score_group_tables(INTEGER(st_par), INTEGER(st_apar), INTEGER(gt_par), INTEGER(gt_apar));
	
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
