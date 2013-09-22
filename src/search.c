/**
 * Triplex search package
 * Search algorithm
 *
 * @author  Matej Lexa, Tomas Martinek, Jiri Hon
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

/* Status variable flags */
#define STAT_NONE     0
#define STAT_QUALITY  1
#define STAT_MINLEN   2
#define STAT_EXPORT   4


double LAMBDA[NUM_TRI_TYPES] = {0.71, 0.71, 0.67, 0.67, 0.71, 0.71, 0.67, 0.67};
double MI[NUM_TRI_TYPES]     = {5.88, 5.88, 6.05, 6.05, 5.88, 5.88, 6.05, 6.05};
const int TAB_STRAND[NUM_TRI_TYPES] = {0, 0, 1, 1, 1, 1, 0, 0};

typedef enum
{// States for diagonal analysis in get_triplex_regions
	S_AD_INIT,
	S_AD_TRIPLEX,
	S_AD_MIN_GAP,
	S_AD_GAP
} ad_states_t;


/** Function prototypes **/
void export_data(t_diag diag, int tri_type, int offset);
void search(
	char *piece, int piece_l, int offset, int n_antidiag, int max_bonus,
	t_diag *diag, t_params *params, t_penalization *pen, prog_t *pb
);

void print_score_array(t_diag* ptr, int size, int border);
void print_rule_array(t_diag* ptr, int size, int border);
void print_status_array(t_diag* ptr, int size, int border);



/**
 * P-value
 * @param score Triplex score
 * @param tri_type Triplex type
 * @return P-value
 */
static inline double p_value(int score, int tri_type)
{
	return 1-exp(-exp(-LAMBDA[tri_type]*(score-MI[tri_type])));
}


/**
 * Search triplex in DNA sequence
 * @param dna encoed DNA sequence, @see encode_bases
 * @param chunk Interval list of chunks divided by N or - symbols
 * @param params Algorithm options
 * @param pen Custom penalizations
 * @param pbw Progress bar width
 */
void main_search(seq_t dna, intv_t *chunk, t_params *params, t_penalization *pen, int pbw)
{
	int chunk_len, npieces, delta, piece_l, last_piece_l, min_score;
	
	t_diag *diag = malloc(3*MAX_PIECE_SIZE * sizeof(t_diag));
	//static t_diag diag[3*MAX_PIECE_SIZE]; // One MAX_PIECE_SIZE extra for piece_overlap
	
	// Maximal bonus per match
	int max_bonus = get_max_bonus(params->tri_type, pen->iso_stay);
	
	min_score = get_min_score(params->p_val, params->tri_type);
	if (min_score > params->min_score)
	// Use minimal score deduced from P-value for better performance
		params->min_score = min_score;
	
	/* NOTE: Maybe now could be all the filtration conditions considering
	 * p_val removed from search() and get_max_score()
	 * as this is now well represented by min_score */
	
	// Number of antidiagonals per triplex
	int n_antidiag = get_n_antidiag(
		max_bonus, pen->insertion, params->max_len, params->min_score,
		params->max_loop
	);
	
	int pieces_overlap = n_antidiag;
	
	/* Initialize progress bar structure */
	prog_t pb = {0, dna.len, pbw};
	
	if (pb.max >= PB_SHOW_LIMIT)
	/* Draw progress bar initially */
		set_txt_progress_bar(&pb, 0);
	
	while (chunk != NULL)
	{
		chunk_len = chunk->end - chunk->start + 1;
		npieces = ceil(chunk_len / (double) MAX_PIECE_SIZE);
		delta = MAX_PIECE_SIZE;
		
		piece_l = MAX_PIECE_SIZE + pieces_overlap;
		last_piece_l = chunk_len - (npieces-1)*MAX_PIECE_SIZE;
		
		/* If last piece is shorter than overlap, then remove it from computation
		 * because previous piece (if exist) calculates its */
		if ((last_piece_l <= pieces_overlap) && (npieces > 1))
		{
			npieces--;
			last_piece_l = chunk_len - (npieces-1)*MAX_PIECE_SIZE;
			/* NOTE: should be same as
			 * last_piece_l = MAX_PIECE_SIZE + last_piece_l */
		}
		
		for (int j = 0, piece_offset = chunk->start;
		     j < npieces;
		     j++, piece_offset += delta)
		{// Iterate through pieces
			
			if (j == npieces-1) piece_l = last_piece_l;
			
			// Diag structure initialization
			for (int i = 0; i < 2*piece_l; i++)
			{
				diag[i].score = 0;
				diag[i].max_score = 0;
				diag[i].bound = 0;
				diag[i].twist = 90;
				diag[i].dtwist = 0;
				diag[i].status = STAT_NONE;
				diag[i].start.diag = i;
				diag[i].start.antidiag = (((params->min_loop+i) % 2) == 0) ? params->min_loop+1 : params->min_loop+2;
				diag[i].max_score_pos.diag = diag[i].start.diag;
				diag[i].max_score_pos.antidiag = diag[i].start.antidiag;
				diag[i].indels = 0;
				diag[i].max_indels = 0;
				diag[i].dp_rule = DP_MISMATCH;
			}  
			search(dna.seq + piece_offset, piece_l, piece_offset, n_antidiag, max_bonus, diag, params, pen, &pb);
		}
		chunk = chunk->next;
	}
	
	free(diag);
	
	if (pb.max >= PB_SHOW_LIMIT)
		Rprintf("\n");
}


/**
 * Deduce minimal score from P-value treshold
 * NOTE P-value must be monotonic decreasing function with respect to score
 * @param pvalue P-value treshold
 * @param type Triplex type
 * @return Minimal score
 */
int get_min_score(double pvalue, int type)
{
	int score = 1;
	
	while (p_value(score, type) > pvalue)
		score++;
	
	return score;
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
	
	for (i=border; i<=size-border; i=i+2)
	{ 
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
	
	for (i = border; i <= size-border; i = i+2)
	{ 
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
	
	for (i=border; i<=size-border; i++)
	{ 
		if (ptr[i].status == STAT_NONE) Rprintf(" ");
		if (ptr[i].status == STAT_EXPORT) Rprintf("*");
		if (ptr[i].status == (STAT_MINLEN|STAT_QUALITY)) Rprintf("|");
		if (ptr[i].status == STAT_QUALITY) Rprintf(".");
	}
	Rprintf("\n");
}


/**
 * Convert diagonal and antidiagonal coordinates
 * into row index starting from bottom
 * NOTE Make sure that d >= ad.
 * @param ad Antidiagonal number
 * @param d  Diagonal number
 * @return Row index starting from bottom
 */
static inline int d_to_start(int ad, int d)
{
/* Example: ad = 2, d = 5 -> result = 1

     0   1   2   3   4   5
   +---+---+---+---+---+---+
5  |   |   |   |   |   |   |
   +---+---+---+---+---+---+
4  |   |   |   |   |   | 10|
   +---+---+---+---+---+---+
3  |   |   |   |   | 8 | 9 |
   +---+---+---+---+---+---+
2  |   |   |   | 6 | 7 |   |
   +---+---+---+---+---+---+
1  |===|===| 4 | 5 |===|===|
   +---+---+---+---+---+---+
0  |   | 2 | 3 |   |   |   |
   +---+---+---+---+---+---+
  O   1   2   3   4   5 
  
*/
  return (d - ad) / 2;
}


/**
 * Convert diagonal and antidiagonal coordinates
 * into column index starting from left
 * NOTE Make sure that d > ad.
 * @param ad Antidiagonal number
 * @param d  Diagonal number
 * @return Column index starting from left
 */
static inline int d_to_end(int ad, int d)
{
/* Example: ad = 2, d = 8 -> result = 4

     0   1   2   3   4   5
   +---+---+---+---+---+---+
5  |   |   |   |   | # |   |
   +---+---+---+---+---+---+
4  |   |   |   |   | # | 10|
   +---+---+---+---+---+---+
3  |   |   |   |   | 8 | 9 |
   +---+---+---+---+---+---+
2  |   |   |   | 6 | 7 |   |
   +---+---+---+---+---+---+
1  |   |   | 4 | 5 | # |   |
   +---+---+---+---+---+---+
0  |   | 2 | 3 |   | # |   |
   +---+---+---+---+---+---+
  O   1   2   3   4   5 
  
*/
	return (d - ad - 1)/2 + ad;
}


/**
 * Create triplex region interval
 * @param start Start diagonal
 * @param end End diagonal
 * @param d_overlap Number of diagonals to overlap
 * @param ad Antidiagonal index
 * @param d_first First diagonal index
 * @param d_last Last diagonal index
 * @return Triplex region interval
 */
static inline intv_t *triplex_region(int start, int end, int d_overlap, int ad, int d_first, int d_last)
{
	start -= d_overlap;
	if (start < d_first)
		start = d_first;
	
	end += d_overlap;
	if (end > d_last)
		end = d_last;
	
	return new_intv(d_to_start(ad, start), d_to_end(ad, end));
}


/**
 * Get intervals whic still need computation
 * @param ad Current antidiagonal index
 * @param n_adiag Number of antidiagonal
 * @param adiag Antidiagonal  
 * @param region Regions to analyze on diagonal
 * @param treshold Minimal score for intervals that still need further computation 
 * @return Intervals which still need computation
 */
intv_t *get_triplex_regions(
	int ad, int n_adiag, t_diag *diag,
	intv_t *region, int treshold)
{
	/* Illustration of diagonal and antidiagonal numbers
	 * 
	 * n_adiag = 6 (number of diagonals)
	 * ad = 2 (current antidiagonal)
	
	     R   O   T   O   R   S
	   +---+---+---+---+---+---+
	S  |   |   |   |   |   |   |
	   +---+---+---+---+---+---+
	R  |   |   |   |   |   | 10|
	   +---+---+---+---+---+---+
	O  |   |   |   |   | 8 | 9 |
	   +---+---+---+---+---+---+
	T  |   |   |   | 6 | 7 |   |
	   +---+---+---+---+---+---+
	O  |   |   | 4 | 5 |   |   |
	   +---+---+---+---+---+---+
	R  |   | 2 | 3 |   |   |   |
	   +---+---+---+---+---+---+
	  O   1   2   3   4   5  <- antidiagonal numbers
	
	 * Diagonal numbers are presented in the table cells
	 */
	
	/* How many diagonals should be added to both sides
	 * of possible triplex interval */
	int d_overlap = n_adiag - ad;
	
	/* Minimal gap between two triplex forming regions */
	int min_gap = 3*d_overlap;
	
	int d, gap_len, d_last, d_first, start, end;
	ad_states_t state;
	intv_t *tmp;
	
	// Create first empty interval as a list header
	intv_t header = {0, 0, NULL};
	intv_t *last = &header;
	
#ifndef NDEBUG
	int triplex = 0;
#endif
	
	while (region != NULL)
	{
		d_first = ad + 2*region->start; // First diagonal
		d_last = 2*(region->end + 1) - ad; // Last diagonal
		
		state = S_AD_INIT;
		start = d_first; // Interval start diagonal
		end = d_last; // Interval end diagonal
		gap_len = 0;
	
		/* Finite-state automata, that looks for regions
		 * with score higher or equal to treshold.
		 * Such regions are extended by d_overlap on both sides if possible.
		 * See triplex_region function for details. */
		for (d = d_first; d <= d_last; d++)
		{
#ifndef NDEBUG
			if (diag[d].score >= treshold)
				triplex++;
#endif
			switch (state)
			{
				case S_AD_INIT:
				// Initial decision
					if (diag[d].score >= treshold)
					{
						start = d;
						state = S_AD_TRIPLEX;
					}
					break;
				case S_AD_TRIPLEX:
				// Triplex forming region
					if (diag[d].score < treshold)
					{
						end = d - 1;
						gap_len = 1;
						state = S_AD_MIN_GAP;
					}
					break;
				case S_AD_MIN_GAP:
				// Check if the gap is at least min_gap long
					if (diag[d].score < treshold)
					{
						gap_len++;
						
						if (gap_len == min_gap)
							state = S_AD_GAP;
					}
					else
						state = S_AD_TRIPLEX;
					break;
				case S_AD_GAP:
				// The gap is long enough
					if (diag[d].score >= treshold)
					{// Export triplex interval
						last->next = triplex_region(start, end, d_overlap, ad, d_first, d_last);
						last = last->next;
						start = d;
						state = S_AD_TRIPLEX;
					}
					break;
			}
		}
		if (state == S_AD_TRIPLEX ||
			state == S_AD_MIN_GAP ||
			state == S_AD_GAP)
		{// Export last triplex interval
			last->next = triplex_region(start, end, d_overlap, ad, d_first, d_last);
			last = last->next;
		}
		tmp = region;
		region = region->next;
		free(tmp);
	}
	
#ifndef NDEBUG
	printf("Possible triplexes: %d\n", triplex);
	
	/* Debug print */
	intv_t *intv = header.next;
	
	int width = 0;
	int count = 0;
	while (intv != NULL)
	{
		width += intv->end - intv->start + 1;
		count++;
		intv = intv->next;
	}
	printf("Adiag: %d, treshold: %d, number of intervals: %d, average interval width: %g\n",
	       ad, treshold, count, (double) width/count);
#endif
	return header.next;
}


/**
 * Search for triplexes in given piece
 * @param piece Given sequence
 * @param piece_l Length of given sequence
 * @param offset Offset from the real start of sequence
 * @param n_antidiag Number of antidiagonals to compute
 * @param max_bonus Maximal bonus per match
 * @param diag t_diag array used to search for triplexes
 * @param params Application parameters
 * @param pen Penalization scores
 * @param pb Progress bar
 */
void search(
	char *piece, int piece_l, int offset, int n_antidiag, int max_bonus,
	t_diag *diag, t_params *params, t_penalization *pen, prog_t *pb)
{
	int i, ad, d, length, treshold, d_count, d_under_tres, ad_start;
	double tres_ratio;
	
	// Starting antidiagonal
	ad_start = params->min_loop + 1;
	
	if (piece_l < n_antidiag)
		n_antidiag = piece_l;
	
	intv_t *triplex_regions = new_intv(0, piece_l - 1);
	intv_t *intv = NULL;
	
	FILE *fd = fopen("dp_matrix_vis_nooptim.dat", "a");
	
	printf("n_antidiag: %d\n", n_antidiag);
	
	/* ad = antidiagonal number */
	for (ad = ad_start; ad < n_antidiag; ad++)
	{
		intv = triplex_regions;
		d_count = 0;
		d_under_tres = 0;
		
		/* Minimal score to still have a chance to satisfy min_score param
		 * at the maximal antidiagonal. */
		treshold = params->min_score - (n_antidiag - ad + 1)/2 * max_bonus;
		
		while (intv != NULL)
		{
			for (i = ad + intv->start, d = ad + 2*intv->start + 1;
			     i <= intv->end;
			     i++, d += 2)
			{
				fprintf(fd, "%d %d\n", i + offset, i-ad + offset);
				
				/* Max score and length calcualtion */
				get_max_score(
					piece[i], piece[i-ad],
					&diag[d-1], &diag[d], &diag[d+1],
					d, ad, params->tri_type, params->max_loop, pen
				);
				length = get_length(
					diag[d].start.antidiag,
					diag[d].max_score_pos.antidiag,
					diag[d].max_indels
				);
				diag[d].status = (length >= params->min_len) ? diag[d].status|STAT_MINLEN : diag[d].status&(~STAT_MINLEN);
				
				/* Actual score satisfies the required quality */
				if (diag[d].score >= params->min_score)
				{
					diag[d].status |= STAT_QUALITY;
					/* If triplex can not continue, then export */
					if ((diag[d].status & STAT_MINLEN) && ((d == (ad+1)) || (d == (2*piece_l-ad-1))))
					{
						diag[d].status = STAT_EXPORT;
						if (p_value(diag[d].max_score, params->tri_type) <= params->p_val)
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
						if (p_value(diag[d].max_score, params->tri_type) <= params->p_val)
						{
							export_data(diag[d], params->tri_type, offset);
						}
						diag[d].max_score = 0;
					}
					else {
						diag[d].status = STAT_NONE;
					}
				}
				d_count++;
				
				if (diag[d].score < treshold)
					d_under_tres++;
			}
			intv = intv->next;
		}
		
		tres_ratio = (double) d_under_tres / d_count;
		
		if (0)//tres_ratio >= TRES_RATIO)
		{
			triplex_regions = get_triplex_regions(ad, n_antidiag, diag, triplex_regions, treshold);
			
// 			intv = tr;
// 			d_in_regions = 0;
// 			
// 			while (intv != NULL)
// 			{
// 				for (i = ad + intv->start; i <= intv->end; i++)
// 					d_in_regions++;
// 				
// 				intv = intv->next;
// 			}
// 			real_ratio = (double) (d_count - d_in_regions) / d_count;
// 			
// 			fprintf(fd, "%d %g %g\n", ad, tres_ratio, real_ratio);
			
#ifndef NDEBUG
			int triplex = 0;
			for (d = ad; d <= 2*piece_l - ad; d++)
			{
				if (diag[d].score >= treshold)
					triplex++;
			}
			printf("Real possible triplexes: %d\n", triplex);
#endif
		}
	}
	// Free last version of triplex regions
	free_intv(triplex_regions);
	
	fclose(fd);
	
	if (pb->max >= PB_SHOW_LIMIT)
	{// Redraw progress bar
// 		ad_perc = (ad+1 - ad_start) / ad_width;
// 		set_txt_progress_bar(pb, offset + ad_perc * piece_l);
		set_txt_progress_bar(pb, offset + piece_l);
	}
	
	/* Print nonexported triplexes */
	for (i = 1; i < (2*piece_l); i++)
	{
		if ((diag[i].status & STAT_QUALITY) && (diag[i].status & STAT_MINLEN))
		{      
			if (p_value(diag[i].max_score, params->tri_type) <= params->p_val)
				export_data(diag[i], params->tri_type, offset);
		}
	}
}
