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

typedef struct
{// Progress bar object
	double min;
	double max;
	int width;
} prog_t;

//prog_t txt_progress_bar(double min, double max, int width);
void set_txt_progress_bar(prog_t *pb, double value);


#endif // PROGRESS_H
