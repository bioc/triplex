/**
 * Triplex package
 * Header file for interval list
 *
 * @author  Jiri Hon
 * @date    2013/03/20
 * @file    interval.h
 * @package triplex
 */

#ifndef INTERVAL_H
#define INTERVAL_H

typedef struct intv
{
	int start;
	int end;
	struct intv *next;
} intv_t;

intv_t *new_intv(int start, int end);
void free_intv(intv_t *intv);
void print_intv(intv_t *intv);

#endif // INTERVAL_H
