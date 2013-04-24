/**
 * Triplex package
 * Progress bar header file
 *
 * @author  Jiri Hon
 * @date    2013/02/13
 * @file    progress.h
 * @package triplex
 */

#ifndef PROGRESS_H
#define PROGRESS_H

#define PB_SHOW_LIMIT 1000000

typedef struct
{// Progress bar object
	double min;
	double max;
	int width;
} prog_t;

void set_txt_progress_bar(prog_t *pb, double value);

#endif // PROGRESS_H
