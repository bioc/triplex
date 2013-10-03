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
MIN_SCORE     = 1
P_VALUE       = 2
MIN_LEN       = 3
MAX_LEN       = 4
MIN_LOOP      = 5
MAX_LOOP      = 6
LAMBDA_PAR_P  = 7
LAMBDA_PAR_E  = 8
LAMBDA_APAR_P = 9
LAMBDA_APAR_E = 10
MU_PAR_P      = 11
MU_PAR_E      = 12
MU_APAR_P     = 13
MU_APAR_E     = 14
RN_PAR_P      = 15
RN_PAR_E      = 16
RN_APAR_P     = 17
RN_APAR_E     = 18
DTWIST_PEN    = 19
INS_PEN       = 20
ISO_PEN       = 21
ISO_BONUS     = 22
MIS_PEN       = 23

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
ST_PROKARYOTIC = 0
ST_EUKARYOTIC  = 1
ST_AUTOMATIC   = 2 # Not supported yet


###
## Score/group table size
##
NBASES = 4

###
## Scoring constants
##
TM = -9
TS = 2
TW = 1

###
## Isomorphic group constants
##
IN = 0
IA = 1
IB = 2
IC = 1
ID = 2
IE = 3

###
## Nucleic acids constants
##
A = 1
C = 2
G = 3
T = 4
