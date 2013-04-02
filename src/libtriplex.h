/**
 * Triplex package
 * Header file for search and align algorithm
 * common functions and structures
 *
 * @author  Matej Lexa, Tomas Martinek
 * @date    2012/10/15
 * @file    libtriplex.h
 * @package triplex
 */

#ifndef LIBTRIPLEX_H
#define LIBTRIPLEX_H

#define NUM_TRI_TYPES   8

/* DP rule constants */
#define DP_MATCH        0
#define DP_MISMATCH     1
#define DP_LEFT         2
#define DP_RIGHT        3
#define DP_STOP         4
#define DP_NONE         5
#define DP_MAIN_ADIAG   6

#define A 0
#define C 1
#define G 2
#define T 3

typedef struct
{// DP matrix position represented as diagonal and antidiagonal
        int     diag;
        int     antidiag;
} t_pos;

typedef struct
{// penalization values
  int dtwist;
  int mismatch;
  int insertion;
  int iso_change;
  int iso_stay;
} t_penalization;

typedef struct
{// main application parameters 
  int tri_type;
  int min_score;
  double p_val;
  int min_len;
  int max_len;
  int min_loop;
  int max_loop;
  int max_chunk_size;
} t_params;

typedef struct
{// Data stored for every diagonal
	t_pos   start;          /* Position of the first match */
	t_pos   max_score_pos;  /* Maximal score position */
	int     score;          /* Actual score */
	int     max_score;      /* Maximal score */
	int     status;         /* Status */
	int     bound;          /* Type of the triplex bound */

	int     twist;          /* angle between C1 atoms */
	int     dtwist;         /* Change in angle between subsequent triplets */
	int     dp_rule;        /* Previous position: 0 - match, 1 - mismatch, 2 - left, 3 - right  */
	int     indels;         /* number of indels */
	int     max_indels;     /* number of indels from start position to max score postition */
} t_diag;


/* Convert ASCII characters to DNA nukleotide (A=0,C=1,G=2,T=3) */
char char2nukl(char ch);

/* Convert DNA nukleotide to ASCII characters */
char nukl2char(char nukl);

/* Calcutates the length of a path in DP matrix */
int get_length(int start_antidiag, int end_antidiag, int start_diag, int end_diag);

/* Application of DP rule */
t_diag get_max_score(
	unsigned char a, unsigned char b, t_diag dl, t_diag d, t_diag dr,
	int diag, int antidiag, int tri_type, int max_loop, t_penalization *p
);

/* Replaces special nucleic acid codes allowed in FASTA format into the most suitable nucleotide */
void handle_special(char *a, char *b, int triplex_type, t_diag d, t_penalization *p);

/*  Returns an array of corresponding symbols for n */
unsigned char *get_meaning(char n);

#endif // LIBTRIPLEX_H
