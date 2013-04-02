###
## Constants for triplex search and align R interface
##
## Author: Jiri Hon
## Date: 2012/10/28
## Package: triplex
##

###
## Parameter positions in parameter array
## for passing from/to C function
##
MIN_SCORE  = 1
P_VALUE    = 2
MIN_LEN    = 3
MAX_LEN    = 4
MIN_LOOP   = 5
MAX_LOOP   = 6
LAMBDA_1   = 7
LAMBDA_2   = 8
MI_1       = 9
MI_2       = 10
DTWIST_PEN = 11
INS_PEN    = 12
ISO_PEN    = 13
ISO_BONUS  = 14
MIS_PEN    = 15

###
## Positions in result list from C
##
T_START    = 1
T_END      = 2
T_SCORE    = 3
T_P_VALUE  = 4
T_INSDEL   = 5
T_TYPE     = 6
T_L_START  = 7
T_L_END    = 8
T_STRAND   = 9
