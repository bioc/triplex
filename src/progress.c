/**
 * Triplex package
 * Progress bar
 *
 * @author  Jiri Hon
 * @date    2013/02/13
 * @file    progress.c
 * @package triplex
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "progress.h"

#define PB_EXTRA_CHARS 10

/**
 * Print progress bar state
 * @param pb Progress bar structure pointer
 * @param value Actual progress bar value between min and max
 */
void set_txt_progress_bar(prog_t *pb, double value)
{
	double intw = pb->max - pb->min;
	double percent = (value - pb->min) / intw;
	int width = pb->width - PB_EXTRA_CHARS;
	
	int nchar = (int) (percent * width);
	int nblank = width - nchar;
	
	Rprintf("\r  |");
	for (int i = 0; i < nchar; i++)
		Rprintf("=");
	for (int i = 0; i < nblank; i++)
		Rprintf(" ");
	Rprintf("|%4.f%%", percent * 100);
}
