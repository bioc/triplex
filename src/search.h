/**
 * Triplex package
 * Header file for search algorithm
 *
 * @author  Matej Lexa, Tomas Martinek, Jiri Hon
 * @date    2012/10/15
 * @file    search.h
 * @package triplex
 */

#ifndef SEARCH_H
#define SEARCH_H

#include "search_interface.h"
#include "libtriplex.h"
#include "interval.h"

#define MAX_PIECE_SIZE (10*1024)

/* Treshold ratio for triplex regions analysis start,
 * deduced empirically */
#define TRES_RATIO 0.93

extern double RN[NUM_TRI_TYPES];
extern double LAMBDA[NUM_TRI_TYPES];
extern double MI[NUM_TRI_TYPES];

void main_search(seq_t dna, intv_t *chunk, t_params *params, t_penalization *pen, int pbw);
int get_min_score(double pvalue, int type, int seq_len);

#endif // SEARCH_H
