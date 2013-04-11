/**
 * Triplex search package
 * Search algorithm
 *
 * @author  Matej Lexa, Tomas Martinek
 * @date    2012/10/15
 * @file    search.c
 * @package triplex
 */

#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "search.h"
#include "search_interface.h"
#include "progress.h"


double LAMBDA[NUM_TRI_TYPES] = {0.71, 0.71, 0.67, 0.67, 0.71, 0.71, 0.67, 0.67};
double MI[NUM_TRI_TYPES]     = {5.88, 5.88, 6.05, 6.05, 5.88, 5.88, 6.05, 6.05};
const int TAB_STRAND[NUM_TRI_TYPES] = {0, 0, 1, 1, 1, 1, 0, 0};

/* Status variable flags */
#define STAT_NONE     0
#define STAT_QUALITY  1
#define STAT_MINLEN   2
#define STAT_EXPORT   4
#define PB_SHOW_LIMIT 1000000

/** Function prototypes **/
void export_data(t_diag diag, int tri_type, int offset);
void search(char *piece, int piece_l, int offset, t_diag *diag, t_params *params, t_penalization *pen, prog_t *pb);

void print_score_array(t_diag* ptr, int size, int border);
void print_rule_array(t_diag* ptr, int size, int border);
void print_status_array(t_diag* ptr, int size, int border);


/**
 * Search triplex in DNA sequence
 * @param dna    DNA sequence
 * @param params Algorithm options
 * @param pen    Custom penalizations
 * @param pbw    Progress bar width
 */
void main_search(seq_t dna, t_params params, t_penalization pen, int pbw)
{
	int i, j; 
	
	int pieces_overlap = params.max_loop + 2*params.max_len + 10;
	
	/* Initialize progress bar structure */
	prog_t pb = {0, dna.len, pbw};
	
	Rprintf("Searching for triplex type %d...\n", params.tri_type);
	
	if (pb.max >= PB_SHOW_LIMIT)
	/* Draw progress bar initially */
		set_txt_progress_bar(&pb, 0);
	
	/** Split sequence into chunks divided by n, -, or special IUPAC symbols **/
	char *chunk = strtok(dna.seq, "n-rmwdvhbsyk");
	
	while (chunk != NULL)
	{ 
		int chunk_l = strlen(chunk);
		
		/** If chunk is too short, skip iteration and get next chunk **/
// 		if (chunk_l <= 2*pieces_overlap && chunk_l != dna.len)
// 		{
// 			chunk = strtok(NULL, "nN-");
// 			continue;
// 		}
		
		int chunk_offset = chunk - dna.seq;

		/* Translation from a, c, g, t to 0, 1, 2, 3 */
		encode_bases(chunk);
	
		char *piece = NULL;
		int pieces = ceil(chunk_l/(double)(params.max_chunk_size));
		int delta = (params.max_chunk_size) - pieces_overlap;
		
		int piece_l = (params.max_chunk_size) + pieces_overlap;
		int last_piece_l = chunk_l % piece_l;
		
		/* Diag structure array alocation */
		t_diag *diag = Calloc(2*piece_l, t_diag);
		
		for (j = 0, piece = chunk; j < pieces; j++, piece += delta)
		{
			if (j == pieces-1) piece_l = last_piece_l;
			
			int piece_offset = piece - chunk;
			/* Diag structure initialization */
			for(i = 0; i < 2*piece_l; i++)
			{
				diag[i].score = 0;
				diag[i].max_score = 0;
				diag[i].bound = 0;
				diag[i].twist = 90;
				diag[i].dtwist = 0;
				diag[i].status = STAT_NONE;
				diag[i].start.diag = i;
				diag[i].start.antidiag = (((params.min_loop+i) % 2) == 0) ? params.min_loop+1 : params.min_loop+2;
				diag[i].max_score_pos.diag = diag[i].start.diag;
				diag[i].max_score_pos.antidiag = diag[i].start.antidiag;
				diag[i].indels = 0;
				diag[i].max_indels = 0;
				diag[i].dp_rule = DP_MISMATCH;
			}  
			
			search(piece, piece_l, chunk_offset+piece_offset, diag, &params, &pen, &pb);
		}
		
		Free(diag);
		
		chunk = strtok(NULL, "n-");
	}
	
	if (pb.max >= PB_SHOW_LIMIT)
		Rprintf("\n");
}


/**
 * P-function
 * @param score    Triplex score
 * @param tri_type Triplex type
 * @return P-function value
 */
double p_function(int score, int tri_type)
{
	return 1-exp(-exp(-LAMBDA[tri_type]*(score-MI[tri_type])));
}


/**
 * E-value
 * @param score    Triplex score
 * @param tri_type Triplex type
 * @return E-value
 */
double e_value(int score, int tri_type)
{
	return p_function(score, tri_type);
}


/**
 * P-value
 * @param score    Triplex score
 * @param tri_type Triplex type
 * @return P-value
 */
double p_value(int score, int tri_type)
{
	return 1-exp(-p_function(score, tri_type));
}


/**
 * Export triplex
 * @param diag
 * @param tri_type
 * @param offset
 */
void export_data(t_diag diag, int tri_type, int offset)
{
	int start_ch, end_ch;
	int start_gap, end_gap;
	
	/* calculation of string positions */
	end_ch = (diag.max_score_pos.diag + diag.max_score_pos.antidiag - 1)/2;
	start_ch = end_ch - diag.max_score_pos.antidiag;
	
	end_gap = (diag.start.diag + diag.start.antidiag - 1)/2;
	start_gap = end_gap - diag.start.antidiag;
	
	save_result(
		offset + start_ch + 1,
		offset + end_ch + 1,
		diag.max_score,
		p_value(diag.max_score, tri_type),
		diag.max_indels,
		tri_type,
		offset + start_gap + 1 + 1,  /* correction to loop start character */
		offset + end_gap + 1 - 1,    /* correction to loop end character */
      TAB_STRAND[tri_type]
	);
}


/**
 * Print score array
 * @param ptr
 * @param size
 * @param border
 */
void print_score_array(t_diag* ptr, int size, int border) 
{
	int i;
	for (i=0; i<border; i++) Rprintf(";");
	
	for (i=border; i<=size-border; i=i+2) { 
		Rprintf("%d", ptr[i].score);
		Rprintf(";;");
	}
	Rprintf("\n");
}


/**
 * Print rule array
 * @param ptr
 * @param size
 * @param border
 */
void print_rule_array(t_diag* ptr, int size, int border) 
{
	int i;
	for (i = 0; i < border; i++) Rprintf(";");
	
	for (i = border; i <= size-border; i = i+2) { 
		if (ptr[i].dp_rule == DP_MATCH) Rprintf("|");
		if (ptr[i].dp_rule == DP_MISMATCH) Rprintf("x");
		if (ptr[i].dp_rule == DP_LEFT) Rprintf("\\");
		if (ptr[i].dp_rule == DP_RIGHT) Rprintf("/");
		Rprintf(";;");
	}
	Rprintf("\n");
}


/**
 * Print status array
 * @param ptr
 * @param size
 * @param border
 */
void print_status_array(t_diag* ptr, int size, int border) 
{
	int i;
	for (i=0; i<border; i++) Rprintf(" ");
	
	for (i=border; i<=size-border; i++) { 
		if (ptr[i].status == STAT_NONE) Rprintf(" ");
		if (ptr[i].status == STAT_EXPORT) Rprintf("*");
		if (ptr[i].status == (STAT_MINLEN|STAT_QUALITY)) Rprintf("|");
		if (ptr[i].status == STAT_QUALITY) Rprintf(".");
	}
	Rprintf("\n");
}


/**
 * Search for triplexes in given piece
 * @param piece   Given sequence
 * @param piece_l Length of given sequence
 * @param offset  Offset from the real start of sequence
 * @param diag    t_diag array used to search for triplexes
 * @param params  Application parameters
 * @param pen     Penalization scores
 * @param pb      Progress bar
 */
void search(char *piece, int piece_l, int offset, t_diag *diag, t_params *params, t_penalization *pen, prog_t *pb)
{
	int i, x, d, restore = 0, length, exp_num = 0;
	int x_start, x_limit;
	double x_width, x_perc;
	
	char *piece_orig = Calloc(piece_l, char);
	
	// Starting diagonal depth
	x_start = params->min_loop + 1;
	// Diagonal depth limit
	x_limit = 2*params->max_len + params->max_loop;
	
	if (piece_l < x_limit)
		x_limit = piece_l;
	
	x_width = x_limit - x_start;
	
	memcpy(piece_orig, piece, piece_l);
	
        /* x-characters (corresponds to antidiagonal number) */
	for (x = x_start; x < x_limit; x++)
	{
		d = x+1;
		
		for(i = x; i < piece_l; i++, d+=2)
		{
			/** FASTA special nucleic acid codes will be replaced by most suitable nucleotide **/
			if(piece[i] > 3 || piece[i-x] > 3) {        
			handle_special(&piece[i], &piece[i-x], params->tri_type, diag[d], pen);
			restore = 1;
			}
			
			/* Max score and length calcualtion */
			diag[d] = get_max_score(piece[i], piece[i-x], diag[d-1], diag[d], diag[d+1], d, x, params->tri_type, params->max_loop, pen);
			length = get_length(diag[d].start.antidiag, diag[d].max_score_pos.antidiag, diag[d].start.diag, diag[d].max_score_pos.diag);
			diag[d].status = (length >= params->min_len) ? diag[d].status|STAT_MINLEN : diag[d].status&(~STAT_MINLEN);
			
			/* Actual score satisfies the required quality */
			if(diag[d].score >= params->min_score)
			{
				diag[d].status |= STAT_QUALITY;
				/* If triplex can not continue, then export */
				if ((diag[d].status & STAT_MINLEN) && ((d == (x+1)) || (d == (2*piece_l-x-1))))
				{
					exp_num++;
					diag[d].status = STAT_EXPORT;
					if(p_value(diag[d].max_score, params->tri_type) < params->p_val)
					{
						export_data(diag[d], params->tri_type, offset);                               
					}          
				}
			}
			/* Actual score does not satisfy the required quality */
			else
			{
			/* If quality requirement was satisfied in previous step, then export */
				if(
				(!(diag[d-1].status & STAT_QUALITY)) && (!(diag[d+1].status & STAT_QUALITY)) &&
				((diag[d].status & STAT_QUALITY)) && ((diag[d].status & STAT_MINLEN)))
				{
					diag[d].status = STAT_EXPORT;
					exp_num++;          
					if(p_value(diag[d].max_score, params->tri_type) < params->p_val)
					{
						export_data(diag[d], params->tri_type, offset);
					}
					diag[d].max_score = 0;
				}
				else {
					diag[d].status = STAT_NONE;
				}
			}
			
		} 
		
		if (pb->max >= PB_SHOW_LIMIT)
		{// Redraw progress bar
			x_perc = (x+1 - x_start) / x_width;
			set_txt_progress_bar(pb, offset + x_perc * piece_l);
		}
		
		//print_rule_array(diag, 2*chunk_l, x+1);
		//print_score_array(diag, 2*chunk_l, x+1);
		
		/** Restore original sequence if any special FASTA nucleic acid codes were replaced **/
		if(restore) {
			restore = 0;
			memcpy(piece, piece_orig, piece_l);
		}   
	}  
	
	/************************* Post-Processing **********************************/  
	/* Print nonexported triplexes */
	int nexp_num = 0;
	for(i = 1; i < (2*piece_l); i++) {
		if ((diag[i].status & STAT_QUALITY) && (diag[i].status & STAT_MINLEN)) {      
			if(p_value(diag[i].max_score, params->tri_type) < params->p_val) {
			export_data(diag[i], params->tri_type, offset);        
			}
			nexp_num++;
		}
	}
	Free(piece_orig);
}
