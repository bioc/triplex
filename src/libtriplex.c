/**
 * Triplex package
 * Search and align algorithm common functions
 *
 * @author  Matej Lexa, Tomas Martinek
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

/* special nucleic acid codes */
unsigned char special[][4] = {
  {A, '|'},           // A
  {C, '|'},           // C
  {G, '|'},           // G
  {T, '|'},           // T
  {A, G, '|'},        // Purine - G/A
  {C, T, '|'},        // Pyrimidine - C/T
  {G, T, '|'},        // Ketone - G/T
  {A, C, '|'},        // Amino group - A/C
  {C, G, '|'},        // Strong interaction - C/G
  {A, T, '|'},        // Weak interaction - A/T
  {C, G, T, '|'},     // Not A
  {A, G, T, '|'},     // Not C
  {A, C, T, '|'},     // Not G
  {A, C, G, '|'},     // Not T
};

/*
Tritype
 0: parallel triplex on first strand
 1: parallel triplex on first strand
 2: parallel triplex on second strand
 3: parallel triplex on second strand
 4: antiparallel triplex on second strand
 5: antiparallel triplex on second strand
 6: antiparallel triplex on first strand
 7: antiparallel triplex on first strand
*/

/* Lookup tables */
const int tab_score[NUM_TRI_TYPES][4][4] = {
 /*
  * Parallel
  */

 /* [H][W] */ /*  T   G   C   A */
 /*   A   */  { {-9, -9, -9, -9},
 /*   C   */    {-9,  2, -9, -9},
 /*   G   */    { 2,  1, -9, -9},
 /*   T   */    {-9,  1,  1,  2} },

 /* [W][H] */ /*  A   C   G   T */
 /*   T   */  { {-9, -9,  2, -9},
 /*   G   */    {-9,  2,  1,  1},
 /*   C   */    {-9, -9, -9,  1},
 /*   A   */    {-9, -9, -9,  2} },

 /* [H][W] */ /*  A   C   G   T */
 /*   T   */  { { 2,  1,  1, -9},
 /*   G   */    {-9, -9,  1,  2},
 /*   C   */    {-9, -9,  2, -9},
 /*   A   */    {-9, -9, -9, -9} },

 /* [W][H] */ /*  T   G   C   A */
 /*   A   */  { { 2, -9, -9, -9},
 /*   C   */    { 1, -9, -9, -9},
 /*   G   */    { 1,  1,  2, -9},
 /*   T   */    {-9,  2, -9, -9} },

 /*
  * Antiparallel
  */

 /* [H][W] */ /*  T   G   C   A */
 /*   T   */  { {-9, -9,  1,  2},
 /*   G   */    {-9,  2, -9, -9},
 /*   C   */    {-9, -9, -9,  1},
 /*   A   */    {-9,  1, -9,  2} },

 /* [W][H] */ /*  T   G   C   A */
 /*   T   */  { {-9, -9, -9, -9},
 /*   G   */    {-9,  2, -9,  1},
 /*   C   */    { 1, -9, -9, -9},
 /*   A   */    { 2, -9,  1,  2} },

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

const int tab_bound[NUM_TRI_TYPES][4][4] = {
 /*
  * Parallel
  */

 /* [H][W] */ /*  T   G   C   A */
 /*   A   */  { { 0,  0,  0,  0},
 /*   C   */    { 0,  1,  0,  0},
 /*   G   */    { 2,  3,  0,  0},
 /*   T   */    { 0,  3,  3,  1} },

 /* [W][H] */ /*  A   C   G   T */
 /*   T   */  { { 0,  0,  2,  0},
 /*   G   */    { 0,  1,  3,  3},
 /*   C   */    { 0,  0,  0,  3},
 /*   A   */    { 0,  0,  0,  1} },

 /* [H][W] */ /*  A   C   G   T */
 /*   T   */  { { 1,  3,  3,  0},
 /*   G   */    { 0,  0,  3,  2},
 /*   C   */    { 0,  0,  1,  0},
 /*   A   */    { 0,  0,  0,  0} },

 /* [W][H] */ /*  T   G   C   A */
 /*   A   */  { { 1,  0,  0,  0},
 /*   C   */    { 3,  0,  0,  0},
 /*   G   */    { 3,  3,  1,  0},
 /*   T   */    { 0,  2,  0,  0} },

 /*
  * Antiparallel
  */

 /* [H][W] */ /*  T   G   C   A */
 /*   T   */  { { 0,  0,  2,  1},
 /*   G   */    { 0,  2,  0,  0},
 /*   C   */    { 0,  0,  0,  1},
 /*   A   */    { 0,  3,  0,  1} },

 /* [W][H] */ /*  T   G   C   A */
 /*   T   */  { { 0,  0,  0,  0},
 /*   G   */    { 0,  2,  0,  3},
 /*   C   */    { 2,  0,  0,  0},
 /*   A   */    { 1,  0,  1,  1} },

 /* [H][W] */ /*  A   C   G   T */
 /*   A   */  { { 1,  0,  3,  0},
 /*   C   */    { 1,  0,  0,  0},
 /*   G   */    { 0,  0,  2,  0},
 /*   T   */    { 1,  2,  0,  0} },

 /* [W][H] */ /*  A   C   G   T */
 /*   A   */  { { 1,  1,  0,  1},
 /*   C   */    { 0,  0,  0,  2},
 /*   G   */    { 3,  0,  2,  0},
 /*   T   */    { 0,  0,  0,  0} },

};

/*
 * Triplet C1-C1-C1 angle
 */
const int tab_twist[NUM_TRI_TYPES][4][4] = {
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


char char2nukl(char ch)
{
	char res;
	ch = tolower(ch);
	switch(ch) {
		case 'r':                   //  purine
		case 'm':                   //  amino group
		case 'w':                   //  weak interaction
		case 'd':                   //  not C
		case 'v':                   //  not T
		case 'h':                   //  not G
		case 'b':                   //  not A
		case 's':                   //  strong interaction
		case 'y':                   //  pyrimidine
		case 'k':                   //  ketone
			res = ch;
			break;
		
		case 'a': res = A; break;   //  adenine
		case 'c': res = C; break;   //  cytosine
		case 'g': res = G; break;   //  guanine
		case 't': res = T; break;   //  thymine    
		
		case 'n':                   //  any nucleotide
		case '-':                   //  insertion of unspecified length
		default:  res = -1;
	}
	return res;
}


/**
 * Get triplex length
 * @param start_antidiag
 * @param end_antidiag
 * @param start_diag
 * @param end_diag
 * @return Triplex length
 */
int get_length(int start_antidiag, int end_antidiag, int start_diag, int end_diag)
{
	int length;
	
	length = (end_antidiag-start_antidiag+2)/2;
	/* diagonal odd to even - increment length */
	if (((start_diag % 2) == 1) && ((end_diag % 2) == 0))
		length++;
	if (length == 0)
		length = 1;
	
	return length;
}


/**
 * Get maximal score
 */
t_diag get_max_score(
	unsigned char a, unsigned char b, t_diag dl, t_diag d, t_diag dr,
	int diag, int antidiag, int tri_type, int max_loop, t_penalization *p)
{
	t_diag res;
	int incscore, mm_score;
	
	/* --------------------------------------------------------------- */
	/* match or mismatch score computation  */
	incscore = tab_score[tri_type][a][b];
		// a trick to have zero or negative scores 
		// for low quality triplets, changed all incscore comparisons to -9
	if (incscore > -9) {  /* match */
		mm_score = d.score + incscore;
		if(d.dp_rule == DP_MATCH) {
			if((tab_bound[tri_type][a][b] != d.bound) && (abs(tab_twist[tri_type][a][b] - d.twist) > p->dtwist ) && (abs(tab_twist[tri_type][a][b] - d.twist + d.dtwist) > p->dtwist)) {      
			mm_score -= p->iso_change;
			}
			else {
			mm_score += p->iso_stay;
			}                                      
		}
	}
	else { /* mismatch */
		mm_score = d.score - p->mismatch;
	}
	
	/* match/mismatch is better */
	if ((mm_score >= dl.score - p->insertion) && (mm_score >= dr.score - p->insertion)) {
		res = d;
		res.score = mm_score;
		res.dp_rule = DP_MISMATCH;
		
		/* match */
		if(incscore > -9) {
			res.dp_rule = DP_MATCH;
			res.bound = tab_bound[tri_type][a][b];
			res.twist = tab_twist[tri_type][a][b];
			res.dtwist = res.twist - d.twist;
			if (mm_score >= res.max_score) {
			res.max_score = mm_score;
			res.max_score_pos.diag = diag;
			res.max_score_pos.antidiag = antidiag;
			res.max_indels = res.indels;
			}
		}
	}
	/* --------------------------------------------------------------- */
	/* insertion or deletion */
	else {
		/* get from left diagonal */
		if (dl.score > dr.score) {
			res = dl;
			res.score = dl.score - p->insertion;
			res.dp_rule = DP_LEFT;
		}
		/* get from right diagonal */
		else {
			res = dr;
			res.score = dr.score - p->insertion;
			res.dp_rule = DP_RIGHT;
		}
		res.indels++;
	}
	/* --------------------------------------------------------------- */
	/* local alignment - only for loop */
	if ((res.score < 0) && (antidiag <= max_loop)) {
		res.score = 0;
		res.max_score = 0;
		res.start.diag = diag;
		res.start.antidiag = antidiag;
		res.max_score_pos.diag = diag;
		res.max_score_pos.antidiag = antidiag;
		res.indels = 0;
		res.max_indels = 0;
	}
	
	/* --------------------------------------------------------------- */
	return res;
}


/**
 * Handle special FASTA symbols
 */
void handle_special(char *a, char *b, int triplex_type, t_diag d, t_penalization *p)
{
	unsigned char *n1 = get_meaning(*a);
	unsigned char *n2 = get_meaning(*b);  
	
	unsigned char i = 0, j = 0, max = - 99, tmp;
	while(n1[i] != '|') {
		j = 0;
		while(n2[j] != '|') {
			tmp = tab_score[triplex_type][n1[i]][n2[j]];
			if((tab_bound[triplex_type][n1[i]][n2[j]] != d.bound) && (abs(tab_twist[triplex_type][n1[i]][n2[j]] - d.twist) > p->dtwist ) && (abs(tab_twist[triplex_type][n1[i]][n2[j]] - d.twist + d.dtwist) > p->dtwist)) {      
				tmp -= p->iso_change;
			}
			else {
				tmp += p->iso_stay;
			}                                            
			if(tmp > max) {
				max = tmp;
				*a = n1[i];
				*b = n2[j];
			}      
			j++;
		}
		i++;
	}  
}


/**
 * Get meaning of special FASTA symbol
 */
unsigned char *get_meaning(char n)
{
	switch(n) {
		case 0:   return special[0];
		case 1:   return special[1];
		case 2:   return special[2];
		case 3:   return special[3];
		case 'r': return special[4];
		case 'y': return special[5];
		case 'k': return special[6];
		case 'm': return special[7];
		case 's': return special[8];
		case 'w': return special[9];
		case 'b': return special[10];
		case 'd': return special[11];
		case 'h': return special[12];
		case 'v': return special[13];
		default:  return special[0];
	}
}
