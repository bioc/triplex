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

#define MAX_CHUNK_SIZE (5*1024*1024)


extern double LAMBDA[NUM_TRI_TYPES];
extern double MI[NUM_TRI_TYPES];

void main_search(seq_t dna, t_params params, t_penalization pen, int pbw);
double p_function(int score, int tri_type);
double p_value(int score, int tri_type, int seq_len);

#endif // SEARCH_H
