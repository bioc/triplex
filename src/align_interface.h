/**
 * Triplex package
 * C interface header file for triplex alignment
 *
 * @author  Jiri Hon
 * @date    2012/10/31
 * @file    align_interface.h
 * @package triplex
 */

#ifndef ALIGN_INTERFACE_H
#define ALIGN_INTERFACE_H

#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>


void Aprintf(char ch);
SEXP triplex_align(SEXP seq, SEXP type, SEXP params);

#endif // ALIGN_INTERFACE_H
