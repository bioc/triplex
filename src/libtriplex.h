/**
 * Triplex package
 * Header file for search and align algorithm
 * common functions and structures
 *
 * @author  Matej Lexa, Tomas Martinek, Jiri Hon
 * @date    2012/10/15
 * @file    libtriplex.h
 * @package triplex
 */

#ifndef LIBTRIPLEX_H
#define LIBTRIPLEX_H

#include <stddef.h>
#include <stdint.h>

#include "interval.h"

// #define NDEBUG

#define NUM_TRI_TYPES   8
#define NUM_SEQ_TYPES   2

/* DP rule constants */
#define DP_MATCH        0
#define DP_MISMATCH     1
#define DP_LEFT         2
#define DP_RIGHT        3
#define DP_STOP         4
#define DP_NONE         5
#define DP_MAIN_ADIAG   6

#define ASCII_LOW 128
#define NBASES 4

#define A 0
#define C 1
#define G 2
#define T 3

#define INVALID_CHAR -1

typedef enum
{// Enumeration for sequence types
	ST_PR = 0,
	ST_EU,
	ST_AU // Not supported yet
} seqtype_t;

typedef struct
{// DP matrix position represented as diagonal and antidiagonal
	int diag;
	int antidiag;
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
} t_params;

typedef struct
{// Data stored for every diagonal
	t_pos start;          /* Position of the first match */
	t_pos max_score_pos;  /* Maximal score position */
	
	uint8_t bound;        /* Type of the triplex bound */
	uint8_t twist;        /* angle between C1 atoms */
	int8_t dtwist;        /* Change in angle between subsequent triplets */
	uint8_t status;       /* Status */
	
	int16_t score;        /* Actual score */
	int16_t max_score;    /* Maximal score */

	uint8_t dp_rule;      /* Previous position: 0 - match, 1 - mismatch, 2 - left, 3 - right  */
	uint8_t indels;       /* number of indels */
	uint8_t max_indels;   /* number of indels from start position to max score postition */
} t_diag;


typedef struct
{// Structure for decoded sequence
	char *seq;
	int len;
	int type;
} seq_t;

extern char CHAR2NUKL[];
extern const char NUKL2CHAR[];
extern int TAB_SCORE[NUM_TRI_TYPES][NBASES][NBASES];
extern int TAB_GROUP[NUM_TRI_TYPES][NBASES][NBASES];
extern const int COMP[NBASES];

void init_CHAR2NUKL_table();

/* Convert ASCII characters to DNA nukleotide (A=0,C=1,G=2,T=3) */
void encode_bases(seq_t dna);

/* Get chunk intervals without N or - symbols */
intv_t *get_chunks(seq_t dna);

int get_max_bonus(int type, int iso_stay_bonus);

/* Get number of antidiagonals to process */
int get_n_antidiag(
	int max_bonus, int ins_pen, int max_len, int min_score,
	int max_loop
);

/* Convert DNA nukleotide to ASCII characters */
char nukl2char(char nukl);

/* Calcutates the length of a path in DP matrix */
int get_length(int start_antidiag, int end_antidiag, int insertions);

/* Application of DP rule */
void get_max_score(
	unsigned char a, unsigned char b, t_diag *dl, t_diag *d, t_diag *dr,
	int diag, int antidiag, int tri_type, int max_loop, t_penalization *p
);

#endif // LIBTRIPLEX_H
