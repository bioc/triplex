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
LAMBDA_1_P = 7
LAMBDA_1_E = 8
LAMBDA_2_P = 9
LAMBDA_2_E = 10
MI_1_P     = 11
MI_1_E     = 12
MI_2_P     = 13
MI_2_E     = 14
RN_1_P     = 15
RN_1_E     = 16
RN_2_P     = 17
RN_2_E     = 18
DTWIST_PEN = 19
INS_PEN    = 20
ISO_PEN    = 21
ISO_BONUS  = 22
MIS_PEN    = 23

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

###
## Prokaryotic vs. eukaryotic
##
S_AUTO = 0
S_PROKARYOTIC = 1
S_EUKARYOTIC = 2
