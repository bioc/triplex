/**
 * Triplex package
 * Shared library initialization
 *
 * @author  Jiri Hon
 * @date    2012/10/15
 * @file    R_init_triplex.c
 * @package triplex
 */

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "search_interface.h"
#include "align_interface.h"
#include "libtriplex.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

// .Call entry points definitions
static const R_CallMethodDef callMethods[] =
{
/* algorithm.c */
	CALLMETHOD_DEF(triplex_search, 9),
/* triplex_align.c */
	CALLMETHOD_DEF(triplex_align, 7),
	{NULL, NULL, 0}
};


/**
 * Shared library initialization function
 * @param info Info about shared library
 */
void R_init_triplex(DllInfo *info)
{
	init_CHAR2NUKL_table();
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	return;
}
