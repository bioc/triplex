/**
 * Triplex package
 * Interval list
 *
 * @author  Jiri Hon
 * @date    2013/03/20
 * @file    interval.c
 * @package triplex
 */

#include <Rinternals.h>
#include <stdlib.h>
#include "interval.h"


/**
 * Create new interval
 * @param start Interval start
 * @param end Interval end
 * @return Pointer to new interval
 */
intv_t *new_intv(int start, int end)
{
	intv_t *intv = malloc(sizeof(intv_t));
	if (intv == NULL)
		error("Failed to allocate memory for new interval.");
	
	intv->start = start;
	intv->end = end;
	intv->next = NULL;
	return intv;
}


/**
 * Free all linked intervals
 * @param intv First interval
 */
void free_intv(intv_t *intv)
{
	intv_t *tmp;
	while (intv != NULL)
	{
		tmp = intv;
		intv = intv->next;
		free(tmp);
	}
}


/**
 * Print all intervals for debugging
 * @param intv First interval
 */
void print_intv(intv_t *intv)
{
	if (intv == NULL)
		Rprintf("No intervals.\n");
	
	while (intv != NULL)
	{
		Rprintf("%d, %d\n", intv->start, intv->end);
		intv = intv->next;
	}
}
