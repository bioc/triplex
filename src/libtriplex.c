/**
 * Triplex package
 * Search and align algorithm common functions
 *
 * @author  Matej Lexa, Tomas Martinek, Jiri Hon
 * @date    2012/10/15
 * @file    libtriplex.c
 * @package triplex
 */

#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "libtriplex.h"


/* Translation table from ascii symbols to
 * internal representation of DNA bases */
char CHAR2NUKL[ASCII_LOW];

/* Translation table from internal representation
 * of DNA bases to corresponding ASCII symbols */
const char NUKL2CHAR[NBASES] =
{
	[A] = 'a', [C] = 'c',
	[G] = 'g', [T] = 't'
};

/* Translation table from internal representation
 * of DNA bases to complement DNA bases */
const int COMP[NBASES] =
{
	[A] = T, [C] = G,
	[G] = C, [T] = A 
};

/* Table indicating which symbols should be cut off
 * from input DNA sequence, especially IUPAC symbols */
const int CHUNKCHAR[ASCII_LOW] =
{
	['n'] = 1, ['-'] = 1, ['r'] = 1, ['m'] = 1,
	['w'] = 1, ['d'] = 1, ['v'] = 1, ['h'] = 1,
	['b'] = 1, ['s'] = 1, ['y'] = 1, ['k'] = 1
};


/* Possible triplex types
---------------------------------------------------

==== Type 0 ====                   ==== Type 1 ====
5' ---------LOOP                   LOOP--------- 3'
    | | | |           Parallel          | | | |
3' --------         first strand        -------- 5'
    : : : :                             : : : :
3' ---------LOOP                   LOOP--------- 5'

==== Type 2 ====                   ==== Type 3 ====
5' ---------LOOP                   LOOP--------- 3'
    : : : :          Parallel           : : : :
5' --------        second strand        -------- 3'
    | | | |                             | | | |
3' ---------LOOP                   LOOP--------- 5'

==== Type 4 ====                   ==== Type 5 ====
5' ---------                            -------- 3'
    | | | |         Antiarallel         | | | |
3' ---------LOOP   second strand   LOOP--------- 5'
    : : : :                             : : : :
5' ---------LOOP                   LOOP--------- 3'

==== Type 6 ====                   ==== Type 7 ====
3' ---------LOOP                   LOOP--------- 5'
    : : : :         Antiparallel        : : : :
5' ---------LOOP    first strand   LOOP--------- 3'
    | | | |                             | | | |
3' ---------                            -------- 5'

---------------------------------------------------
*/


/* Dynamic algorithm principle
 * is based on palindrome search

Example
Scoring: +2 match, -7 mismatch, -9 ins
Max loop: 1

     R   O   T   O   R   S
   +---+---+---+---+---+---+
S  |   |   |   |   | O | O |
   +---+---+---+---+---+---+
R  |   |   |   | O | 0 | 0 |
   +---+---+---+---+---+---+
O  |   |   | O | 0 | 0 |-7 |
   +---+---+---+---+---+---+
T  |   | 0 | 0 | 0 |-7 |-7 |
   +---+---+---+---+---+---+
O  | 0 | 0 | 0 | 2 |-7 |-14|
   +---+---+---+---+---+---+
R  | 0 | 0 |-7 |-7 | 4 |-5 |
   +---+---+---+---+---+---+
  O   1   2   3   4   5  <- antidiagonal numbers

Result palindrome:

R---O--
|   |  T
R---O--
*/

/* Triplex detection scheme
 * (pseudo-palindromatic structure)

                      --------------
                    /               \
                   |                |
           5' -----Y----------------X-----   3'

Order:           second           first
DP Matrix:        row             column
Score Matrix:    column            row
*/


// Strong triplet score
#define TS 2
// Alternative (weak) triplet score
#define TW 1
// Triplet missmatch indication
#define TM -9

int TAB_SCORE[NUM_TRI_TYPES][NBASES][NBASES] =
{// Tabulated triplet score
 /*
  * Parallel
  */

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {TM, TM, TM, TM},
 /*   C   */    {TM, TS, TM, TM},
 /*   G   */    {TS, TW, TM, TM},
 /*   T   */    {TM, TW, TW, TS} },

 /* [W][H] */ /*  A   C   G   T */
 /*   A   */  { {TM, TM, TS, TM},
 /*   C   */    {TM, TS, TW, TW},
 /*   G   */    {TM, TM, TM, TW},
 /*   T   */    {TM, TM, TM, TS} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {TS, TW, TW, TM},
 /*   C   */    {TM, TM, TW, TS},
 /*   G   */    {TM, TM, TS, TM},
 /*   T   */    {TM, TM, TM, TM} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {TS, TM, TM, TM},
 /*   C   */    {TW, TM, TM, TM},
 /*   G   */    {TW, TW, TS, TM},
 /*   T   */    {TM, TS, TM, TM} },

 /*
  * Antiparallel
  */

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {TM, TM, TW, TS},
 /*   C   */    {TM, TS, TM, TM},
 /*   G   */    {TM, TM, TM, TW},
 /*   T   */    {TM, TW, TM, TS} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {TM, TM, TM, TM},
 /*   C   */    {TM, TS, TM, TW},
 /*   G   */    {TW, TM, TM, TM},
 /*   T   */    {TS, TM, TW, TS} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {TS, TM, TW, TM},
 /*   C   */    {TW, TM, TM, TM},
 /*   G   */    {TM, TM, TS, TM},
 /*   T   */    {TS, TW, TM, TM} },

 /* [W][H] */ /*  A   C   G   T */
 /*   A   */  { {TS, TW, TM, TS},
 /*   C   */    {TM, TM, TM, TW},
 /*   G   */    {TW, TM, TS, TM},
 /*   T   */    {TM, TM, TM, TM} },

};


/* Isomorphic groups named according to triplet table
 * published in Lexa et al. 2011 */
#define IN 0
// Parallel triplets
#define IA 1
#define IB 2
// Antiparallel triplets
#define IC 3
#define ID 4
#define IE 5

int TAB_GROUP[NUM_TRI_TYPES][NBASES][NBASES] =
{// Tabulated isomorphic groups
 /*
  * Parallel
  */

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {IN, IN, IN, IN},
 /*   C   */    {IN, IA, IN, IN},
 /*   G   */    {IB, IB, IN, IN},
 /*   T   */    {IN, IA, IB, IA} },

 /* [W][H] */ /*  A   C   G   T */
 /*   A   */  { {IN, IN, IB, IN},
 /*   C   */    {IN, IA, IB, IA},
 /*   G   */    {IN, IN, IN, IB},
 /*   T   */    {IN, IN, IN, IA} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {IA, IB, IA, IN},
 /*   C   */    {IN, IN, IB, IB},
 /*   G   */    {IN, IN, IA, IN},
 /*   T   */    {IN, IN, IN, IN} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {IA, IN, IN, IN},
 /*   C   */    {IB, IN, IN, IN},
 /*   G   */    {IA, IB, IA, IN},
 /*   T   */    {IN, IB, IN, IN} },

 /*
  * Antiparallel
  */

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {IN, IN, IE, IC},
 /*   C   */    {IN, IE, IN, IN},
 /*   G   */    {IN, IN, IN, ID},
 /*   T   */    {IN, ID, IN, IC} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {IN, IN, IN, IN},
 /*   C   */    {IN, IE, IN, ID},
 /*   G   */    {IE, IN, IN, IN},
 /*   T   */    {IC, IN, ID, IC} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { {IC, IN, ID, IN},
 /*   C   */    {ID, IN, IN, IN},
 /*   G   */    {IN, IN, IE, IN},
 /*   T   */    {IC, IE, IN, IN} },

 /* [W][H] */ /*  A   C   G   T */
 /*   A   */  { {IC, ID, IN, IC},
 /*   C   */    {IN, IN, IN, IE},
 /*   G   */    {ID, IN, IE, IN},
 /*   T   */    {IN, IN, IN, IN} },

};


const int TAB_TWIST[NUM_TRI_TYPES][NBASES][NBASES] =
{// Triplet C1-C1-C1 angle
 /*
  * Parallel
  */

 /* [H][W] */ /*    T    G   C    A */
 /* A */      { {   0,   0,  0,   0},
 /* C */        {   0, 109,  0,   0},
 /* G */        { 126,  75,  0,   0},
 /* T */        {   0,  78, 71, 104} },

 /* [W][H] */ /*   A    C    G    T */
 /*   T   */  { {  0,   0, 126,   0},
 /*   G   */    {  0, 109,  75,  78},
 /*   C   */    {  0,   0,   0,  71},
 /*   A   */    {  0,   0,   0, 104} },

 /* [H][W] */ /*    A   C    G    T */
 /*   T   */  { { 104, 71,  78,   0},
 /*   G   */    {   0,  0,  75, 126},
 /*   C   */    {   0,  0, 109,   0},
 /*   A   */    {   0,  0,   0,   0} },

 /* [W][H] */ /*    T     G    C  A */
 /*   A   */  { { 104,    0,   0, 0},
 /*   C   */    {  71,    0,   0, 0},
 /*   G   */    {  78,   75, 109, 0},
 /*   T   */    {   0,  126,   0, 0} },

 /*
  * Antiparallel
  */

 /* [H][W] */ /*  T     G   C    A */
 /*   T   */  { { 0,    0, 94,  72},
 /*   G   */    { 0,   94,  0,   0},
 /*   C   */    { 0,    0,  0,  72},
 /*   A   */    { 0,  126,  0,  77} },

 /* [W][H] */ /*   T   G   C     A */
 /*   T   */  { {  0,  0,  0,    0},
 /*   G   */    {  0, 94,  0,  126},
 /*   C   */    { 94,  0,  0,    0},
 /*   A   */    { 72,  0, 72,   77} },

 /* [H][W] */ /*  A    C     G   T */
 /*   A   */  { { 77,  0,  126,  0},
 /*   C   */    { 72,  0,    0,  0},
 /*   G   */    {  0,  0,   94,  0},
 /*   T   */    { 72, 94,    0,  0} },

 /* [W][H] */ /*    A    C   G    T */
 /*   A   */  { {  77,  72,  0,  72},
 /*   C   */    {   0,   0,  0,  94},
 /*   G   */    { 126,   0, 94,   0},
 /*   T   */    {   0,   0,  0,   0} },

};


typedef enum
{// States for chunk analysis
	S_CHUNK_INIT,
	S_CHUNK_IN,
	S_CHUNK_OUT
} ch_state_t;


/**
 * Get maximal bonus per diagonal step
 * @param type Triplex type
 * @param iso_stay_bonus Bonus per isogroup stay
 * @return Maximal bonus
 */
int get_max_bonus(int type, int iso_stay_bonus)
{
	int max = 0;
	for (int i = 0; i < NBASES; i++)
	{
		for (int j = 0; j < NBASES; j++)
		{
			if (TAB_SCORE[type][i][j] > max)
				max = TAB_SCORE[type][i][j];
		}
	}
	return max + iso_stay_bonus;
}


/**
 * Get number of antidiagonals to process
 * @param max_bonus Maximal bonus per match
 * @param ins_pen Penalization per insertion 
 * @param max_len Maximal triplex body (stem) length
 * @param min_score Minimal score
 * @param max_loop Maximal loop length
 * @return Number of diagonals to process
 */
int get_n_antidiag(
	int max_bonus, int ins_pen, int max_len, int min_score, int max_loop)
{
	/* Maximal bonus for triplex body.
	 * That means triplex stem length is exactly equal to max_len
	 * and there is maximal amount of insertions to still satisfy
	 * min_score limit */
	int total_bonus = max_bonus * max_len;
	
	// Get score surplus which could absorb penalization for insertions
	int score_surplus = total_bonus - min_score;
	
	int n_ins = 0;
	
	if (score_surplus > 0)
	{// How many insertions can be accepted?
		n_ins = score_surplus / ins_pen;
	}
	// For every match, two antidiagonals must be processed
	return max_loop + 2*max_len + n_ins;
}


char nukl2char(char nukl)
{
	char res;
	switch(nukl){
		case A: res = 'a'; break;
		case C: res = 'c'; break;
		case G: res = 'g'; break;
		case T: res = 't'; break;
		default: res = nukl;
	}
	return res;
}


/**
 * Initialize global CHAR2NUKL translation table
 * Called from @see R_init_triplex
 */
void init_CHAR2NUKL_table()
{
	for (int i = 0; i < ASCII_LOW; i++)
		CHAR2NUKL[i] = INVALID_CHAR;
	
	CHAR2NUKL['a'] = A;
	CHAR2NUKL['c'] = C;
	CHAR2NUKL['g'] = G;
	CHAR2NUKL['t'] = T;
	CHAR2NUKL['n'] = 'n';
	CHAR2NUKL['-'] = '-';
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
}


/**
 * Encode bases to enum values A,C,G,T
 * @param seq DNA sequence
 */
#if 0
void encode_bases(seq_t dna)
{
	char ch;
	for (int i = 0; i < dna.len; i++)
	{
		ch = CHAR2NUKL[tolower(dna.seq[i])];
		if (ch == INVALID_CHAR)
			error("Unsupported symbol '%c' in input sequence.", dna.seq[i]);
		else
			dna.seq[i] = ch;
	}
}
#endif

/**
 * Get chunk intervals without N or - symbols
 * NOTE This function does well only on encoded DNA sequence
 * @param dna encoded DNA sequence @see encode_bases
 * @return Intervals without N or - symbols
 */
intv_t *get_chunks(seq_t dna)
{
	int i, offset = 0;
	ch_state_t state = S_CHUNK_INIT;
	
	// Create first interval as a list header
	intv_t header = {0, 0, NULL};
	intv_t *last = &header;
	
	for (i = 0; i < dna.len; i++)
	{
		switch (state)
		{
			case S_CHUNK_INIT:
			// Initial decision
				if (CHUNKCHAR[(int) dna.seq[i]])
					state = S_CHUNK_OUT;
				else
					state = S_CHUNK_IN;
				break;
			case S_CHUNK_IN:
			// Chunk body
				if (CHUNKCHAR[(int) dna.seq[i]])
				{// Export chunk
					last->next = new_intv(offset, i-1);
					last = last->next;
					state = S_CHUNK_OUT;
				}
				break;
			case S_CHUNK_OUT:
			// Gap between chunks (N or - symbols)
				if (!CHUNKCHAR[(int) dna.seq[i]])
				{// Move offset
					offset = i;
					state = S_CHUNK_IN;
				}
				break;
		}
	}
	if (state == S_CHUNK_IN)
	{// Export last chunk
		last->next = new_intv(offset, i-1);
		last = last->next;
	}
	return header.next;
}


/**
 * Get triplex length
 * @param start_antidiag Starting antidiagonal
 * @param end_antidiag Ending antidiagonal number
 * @param insertions Number of insertions
 * @return Triplex length
 */
int get_length(int start_antidiag, int end_antidiag, int insertions)
{
	return (end_antidiag - start_antidiag - insertions)/2 + 1;
}


/**
 * Get maximal score
 * @param a Symbol from normal sequence
 * @param b Symbol from reversed sequence
 * @param dl Left matrix cell
 * @param d  Current matrix cell
 * @param dr Right matrix cell
 * @param diag Diagonal number
 * @param antidiag Antidiagonal number
 * @param tri_type Triplex type
 * @param max_loop Maximal loop length
 * @param p Penalizations
 */
void get_max_score(
	unsigned char a, unsigned char b, t_diag *dl, t_diag *d, t_diag *dr,
	int diag, int antidiag, int tri_type, int max_loop, t_penalization *p)
{
// 	t_diag res;
	int incscore, mm_score;
	
	/* --------------------------------------------------------------- */
	/* match or mismatch score computation  */
	incscore = TAB_SCORE[tri_type][a][b];
	// a trick to have zero or negative scores 
	// for low quality triplets, changed all incscore comparisons to TM
	if (incscore > TM)
	{// Match
		mm_score = d->score + incscore;
		if (d->dp_rule == DP_MATCH)
		{// Check isomorphic group
			if ((TAB_GROUP[tri_type][a][b] != d->bound) &&
			    (abs(TAB_TWIST[tri_type][a][b] - d->twist) > p->dtwist ) &&
			    (abs(TAB_TWIST[tri_type][a][b] - d->twist + d->dtwist) > p->dtwist))
			{
				mm_score -= p->iso_change;
			}
			else
			{
				mm_score += p->iso_stay;
			}
		}
	}
	else
	{// Mismatch
		mm_score = d->score - p->mismatch;
	}
	
	if ((mm_score >= dl->score - p->insertion) &&
	    (mm_score >= dr->score - p->insertion))
	{// Match/mismatch is better
// 		res = d;
// 		res.score = mm_score;
// 		res.dp_rule = DP_MISMATCH;
		d->score = mm_score;
		d->dp_rule = DP_MISMATCH;
		
		if (incscore > TM)
		{// Match
// 			res.dp_rule = DP_MATCH;
// 			res.bound = TAB_GROUP[tri_type][a][b];
// 			res.twist = TAB_TWIST[tri_type][a][b];
// 			res.dtwist = res.twist - d->twist;
			d->dp_rule = DP_MATCH;
			d->bound = TAB_GROUP[tri_type][a][b];
			d->dtwist = TAB_TWIST[tri_type][a][b] - d->twist;
			d->twist = TAB_TWIST[tri_type][a][b];
			
// 			if (mm_score >= res.max_score)
			if (mm_score >= d->max_score)
			{
// 				res.max_score = mm_score;
// 				res.max_score_pos.diag = diag;
// 				res.max_score_pos.antidiag = antidiag;
// 				res.max_indels = res.indels;
				d->max_score = mm_score;
				d->max_score_pos.diag = diag;
				d->max_score_pos.antidiag = antidiag;
				d->max_indels = d->indels;
			}
		}
	}
	else
	{// Insertion or deletion
		if (dl->score > dr->score)
		{// Get from left diagonal
// 			res = dl;
// 			res.score = dl->score - p->insertion;
// 			res.dp_rule = DP_LEFT;
			*d = *dl;
			d->score = dl->score - p->insertion;
			d->dp_rule = DP_LEFT;
		}
		else
		{// Get from right diagonal
// 			res = dr;
// 			res.score = dr->score - p->insertion;
// 			res.dp_rule = DP_RIGHT;
			*d = *dr;
			d->score = dr->score - p->insertion;
			d->dp_rule = DP_RIGHT;
		}
// 		res.indels++;
		d->indels++;
	}
	
// 	if ((res.score < 0) && (antidiag <= max_loop))
	if ((d->score < 0) && (antidiag <= max_loop))
	{// Local alignment only for loop
// 		res.score = 0;
// 		res.max_score = 0;
// 		res.start.diag = diag;
// 		res.start.antidiag = antidiag;
// 		res.max_score_pos.diag = diag;
// 		res.max_score_pos.antidiag = antidiag;
// 		res.indels = 0;
// 		res.max_indels = 0;
		d->score = 0;
		d->max_score = 0;
		d->start.diag = diag;
		d->start.antidiag = antidiag;
		d->max_score_pos.diag = diag;
		d->max_score_pos.antidiag = antidiag;
		d->indels = 0;
		d->max_indels = 0;
	}
	
// 	return res;
}
