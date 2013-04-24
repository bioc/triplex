/**
 * Triplex package
 * Align algorithm
 *
 * @author  Matej Lexa, Tomas Martinek, Kamil Rajdl, Jiri Hon
 * @date    2012/10/31
 * @file    align.c
 * @package triplex
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "libtriplex.h"

#include "search.h"
#include "align.h"
#include "align_interface.h"

/* Status variable flags */
#define STAT_NONE 0


void search_align(char *piece, int piece_l, t_diag *diag, t_params *params, t_penalization *pen, t_diag** dp_matrix);

t_diag** alloc_matrix(int size);
void free_matrix(t_diag** mat, int size);
void print_matrix(t_diag** mat, char *seq, int l, int compact);
void init_matrix(t_diag** mat, int l);

void print_triplex(t_diag** mat, char* seq, int l);

char rule(int r);


/**
 * Align triplex
 * @param dna    Triplex sequence
 * @param params Algorithm options
 * @param pen    Custom penalizations
 */
void main_align(seq_t dna, t_params params, t_penalization pen)
{      
	int i, verbose_flag = 0;
	
	/* Translation from a, c, g, t to 0, 1, 2, 3 */     
	encode_bases(dna);
	
	/* Diag structure array alocation */
	t_diag *diag = Calloc(2*dna.len, t_diag);
	
	/* Diag structure initialization */
	for(i = 0; i < 2*dna.len; i++)
	{
		diag[i].score = 0;
		diag[i].max_score = 0;
		diag[i].bound = 0;
		diag[i].twist = 90;
		diag[i].dtwist = 0;
		diag[i].status = STAT_NONE;
		diag[i].start.diag = i;
		diag[i].start.antidiag = (((params.min_loop+i) % 2) == 0) ? params.min_loop+1 : params.min_loop+2;
		diag[i].max_score_pos.diag = diag[i].start.diag;
		diag[i].max_score_pos.antidiag = diag[i].start.antidiag;
		diag[i].indels = 0;
		diag[i].max_indels = 0;
		diag[i].dp_rule = DP_MISMATCH;
	}  

	t_diag** mat = alloc_matrix(dna.len);
	
	init_matrix(mat, dna.len);
	
	search_align(dna.seq, dna.len, diag, &params, &pen, mat);
	
	if (verbose_flag) {
		print_matrix(mat, dna.seq, dna.len, 0);
	}
	print_triplex(mat, dna.seq, dna.len);
		
	free_matrix(mat, dna.len);
	Free(diag);
}


/**
 * Prints triplex from the matrix by moving from the bottom right corner
 * of the matrix in a direction given by the rule on the current position until 
 * it reaches main antidiagonals or loop start.
 * @param mat matrix
 * @param seq sequence
 * @param seq_l sequence length
 */
void print_triplex(t_diag** mat, char* seq, int seq_l)
{
	char *body1 = Calloc(seq_l, char);
	char *body2 = Calloc(seq_l, char);
	
	memset(body1, '\0', seq_l);
	memset(body2, '\0', seq_l);
	
	/** Backwards walk from the bottom right corner of the matrix **/
	int i1 = 0, i2 = 0, x, body1_l = 0, body2_l = 0;
	int row = seq_l-1;
	int col = row;
	
	while (mat[row][col].dp_rule != DP_MAIN_ADIAG && mat[row][col].dp_rule != DP_STOP) 
	{
		switch (mat[row][col].dp_rule)
		{
			case DP_MATCH:
				body1[i1++] = nukl2char(seq[seq_l-1-row]);
				body2[i2++] = nukl2char(seq[col]);
				body1_l++;
				body2_l++;
				row--;
				col--;
				break;
			case DP_MISMATCH:
				body1[i1++] = toupper(nukl2char(seq[seq_l-1-row]));
				body2[i2++] = toupper(nukl2char(seq[col]));
				body1_l++;
				body2_l++;
				row--;
				col--;
				break;
			case DP_LEFT:
				body1[i1++] = '-';
				body2[i2++] = nukl2char(seq[col]);
				body2_l++;
				col--;
				break;
			case DP_RIGHT:
				body1[i1++] = nukl2char(seq[seq_l-1-row]);
				body2[i2++] = '-';
				body1_l++;
				row--;
				break;
			default:
				break;
		}
	}
	
	for (x = 0; x < i1; x++)
		Aprintf(body1[x]);
	
	Aprintf('=');
	
	for (x = body1_l; x < (seq_l - body2_l); x++)
		Aprintf(nukl2char(seq[x]));
	
	Aprintf('=');
	
	for (x = i2 - 1; x >= 0; x--)
		Aprintf(body2[x]);
	
	Free(body1);
	Free(body2);
}


/**
 * Maps DP_RULES to chars
 * @param r DP_RULE
 * @return char representing given rule
 */
char rule(int r)
{
	switch(r)
	{
		case DP_MATCH: return '\\';
		case DP_MISMATCH: return 'X';
		case DP_LEFT: return '-';
		case DP_RIGHT: return '|';
		case DP_STOP: return 'S';
		case DP_MAIN_ADIAG: return '/';
		default: return '?';
	}
}


/**
 * Prints DP matrix in CSV format
 * @param mat 2D array representing matrix
 * @param seq sequence
 * @param l sequence length
 * @param compact mode
 */
void print_matrix(t_diag** mat, char* seq, int l, int compact)
{
	if (mat != NULL && seq != NULL)
	{
		int r, c;
		
		if (compact) {
			for(c = 0; c < l; c++) Rprintf(";%c", nukl2char(seq[c]));
			Rprintf("\n");

			for(r = 0; r < l; r++) {
				Rprintf("%c;", nukl2char(seq[l-1-r]));
				for(c = 0; c < l; c++) {
					Rprintf("%c;", rule(mat[r][c].dp_rule));
					//Rprintf("%d;", mat[r][c].max_score);
				}
				Rprintf("\n");
			}
		}
		else {
			Rprintf(" ;");
			for(c = 0; c < l; c++) Rprintf("  %c;", nukl2char(seq[c]));
			Rprintf("\n");

			for(r = 0; r < l; r++) {
				Rprintf("%c;", nukl2char(seq[l-1-r]));
				for(c = 0; c < l; c++) {
					//Rprintf("%d,%c;", mat[r][c].score, rule(mat[r][c].dp_rule));
					Rprintf("  %c;", rule(mat[r][c].dp_rule));
				}
				Rprintf("\n");
				Rprintf("  ");
				for(c = 0; c < l; c++) {
					Rprintf("%3d;", mat[r][c].score);
				}
				Rprintf("\n");
				Rprintf("  ");
				for(c = 0; c < l; c++) {
					Rprintf("%3d;", mat[r][c].start.antidiag);
				}
				Rprintf("\n");
			}
			Rprintf("Score: %d\n", mat[l-1][l-1].score);
			Rprintf("Max indels: %d\n", mat[l-1][l-1].max_indels);
			Rprintf("Indels: %d\n", mat[l-1][l-1].indels);
		}    
	}
}


/**
 * Allocates square 2D array
 * @param size 
 * @return Returns pointer to the allocated array
 */
t_diag** alloc_matrix(int size)
{
	t_diag** mat = (t_diag**)malloc(size*sizeof(t_diag*));  
	if (mat == NULL) {
		error("Not enough space for mat[].");
		return NULL;
	}
	
	int x, y;
	for (x = 0; x < size; x++) {    
		mat[x] = (t_diag*)malloc(size*sizeof(t_diag));      
		if (mat[x] == NULL) {
			for (y = 0; y < x; y++) free(mat[x]);
			error("Not enough space for mat[%d][].", y);
			return NULL;
		}
	}
	return mat;
}


/**
 * Frees given matrix from memory
 * @param mat matrix
 * @param size number of matrixes
 */
void free_matrix(t_diag** mat, int size)
{
	if (mat != NULL)
	{
		for (int i = 0; i < size; i++)
			free(mat[i]);
		
		free(mat);
	}
}


/**
 * Initialize given matrix. 
 * @param mat matrix
 * @param l matrix size
 */
void init_matrix(t_diag** mat, int l)
{
	if (mat != NULL)
	{
		int r, c;
		for (r = 0; r < l; r++)
		{
			for(c = 0; c < l; c++)
			{
				mat[r][c].score = 0;
				mat[r][c].max_score = 0;
				mat[r][c].bound = 0;
				mat[r][c].twist = 90;
				mat[r][c].dtwist = 0;
				mat[r][c].status = STAT_NONE;
				mat[r][c].start.diag = 0;
				mat[r][c].start.antidiag = 0;
				mat[r][c].max_score_pos.diag = mat[r][c].start.diag;
				mat[r][c].max_score_pos.antidiag = mat[r][c].start.antidiag;
				mat[r][c].indels = 0;
				mat[r][c].max_indels = 0;
				/** Two main antidiagonals are initialized differently **/
				if (r == l-1-c || r == l-2-c) {
					mat[r][c].dp_rule = DP_MAIN_ADIAG;
				}        
				else {
					mat[r][c].dp_rule = DP_STOP;
				}
			}  
		}
	}
}


/**
 * Search for triplexes in given piece
 * @param piece Given sequence
 * @param piece_l Length of given sequence
 * @param diag t_diag array used to search for triplexes
 * @param params Application parameters
 * @param pen Penalization scores
 */
void search_align(char *piece, int piece_l, t_diag *diag, t_params *params, t_penalization *pen, t_diag** dp_matrix)
{
	int i, x, d;
	
	for (x = params->min_loop+1; x < piece_l; x++)
	{
		d = x+1;
		for (i = x; i < piece_l; i++, d+=2)
		{
			/* Max score and length calcualtion */
			get_max_score(
				piece[i], piece[i-x],
				&diag[d-1], &diag[d], &diag[d+1],
				d, x, params->tri_type, params->max_loop, pen
			);
			
			/** Resulting dp matrix **/      
			dp_matrix[piece_l-1-(i-x)][i] = diag[d];
		}
	}
}
