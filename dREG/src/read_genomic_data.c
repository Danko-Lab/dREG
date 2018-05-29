/*
 * Pulls genomic data from a bigwig file(s).
 */
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>
#include <assert.h>
#include "read_genomic_data.h"
#include <stdbool.h>
#include <string.h>


#define MERGE_RANGE 1024*1024*20


#define BOOLEAN_ELT(x,__i__)	LOGICAL(x)[__i__]

//#include "bigwiglib.h"

/*
 * set_zoom_params --> Creates a new zoom_params_t data type.
 */
zoom_params_t set_zoom_params(int n_sizes, int* window_sizes, int* half_n_windows) {
  zoom_params_t zoom;
  zoom.window_sizes = window_sizes;
  zoom.half_n_windows = half_n_windows;

  // This is JUST for translating understandable object into a nice vector.
  zoom.n_prev_bins = (int*)R_alloc(n_sizes+1, sizeof(int));
  zoom.n_prev_bins[0]= 0;
  for(int i=1;i<n_sizes;i++) {
    zoom.n_prev_bins[i] = zoom.n_prev_bins[i-1]+2*half_n_windows[i-1];
  }
  zoom.n_prev_bins[n_sizes] = zoom.n_prev_bins[n_sizes-1]+2*half_n_windows[n_sizes-1];

  return(zoom);
}

/*
 * init_genomic_data_point --> Initializes genomic_data_point_t to 0.
 */
void init_genomic_data_point(genomic_data_point_t dp, zoom_params_t zoom) {
  for(int i=0;i<zoom.n_sizes;i++) {
    for(int j=0;j<(2*zoom.half_n_windows[i]);j++) {
      dp.forward[i][j]= 0.0;
      dp.reverse[i][j]= 0.0;
    }
  }
}

/*
 * alloc_genomic_data_point --> Allocates a new genomic_data_point_t.
 */
genomic_data_point_t alloc_genomic_data_point(zoom_params_t zoom) {
  genomic_data_point_t dp;
  dp.forward = (double**)Calloc(zoom.n_sizes,double*);
  dp.reverse = (double**)Calloc(zoom.n_sizes,double*);

  for(int i=0;i<zoom.n_sizes;i++) {
    dp.forward[i] = (double*)Calloc((2*zoom.half_n_windows[i]),double);
    dp.reverse[i] = (double*)Calloc((2*zoom.half_n_windows[i]),double);
  }

  return(dp);
}

/*
 * free_genomic_data_point --> frees a used genomic_data_point_t.
 */
void free_genomic_data_point(genomic_data_point_t dp, zoom_params_t zoom) {
  for(int i=0;i<zoom.n_sizes;i++) {
    Free(dp.forward[i]);
    Free(dp.reverse[i]);
  }

  Free(dp.forward);
  Free(dp.reverse);
}

/*
 * max_dist_from_center --> Returns maximum bounds around a center position.
 */
int max_dist_from_center(int n_sizes, int *window_sizes, int *half_n_windows) {
  int max_bounds=0;
  for(int i=0;i<n_sizes;i++) {
    int curr_bounds= window_sizes[i]*half_n_windows[i];
	if(max_bounds < curr_bounds) max_bounds= curr_bounds;
  }
  return(max_bounds);
}

/*
 * get_bin_number --> Returns the bin that a particular position falls into.
 *
 * Note than since the center base is not included, it needs to be subtracted
 * from the position of interest, if that position falls past the center.
 */
int get_bin_number(int center, int position, int window_size, int half_n_windows) {
  // How to get the bin number with these variables?!
  int left= center-window_size*half_n_windows;
  int right= center+window_size*half_n_windows;
  if(position < left || position > right || position == center) return(-1);
  int dist_from_start= (position<center)?(position-left):(position-left-1); // Because center position isn't included ... have to subtract 1 for windows right of center.
  return((int)floor(((double)dist_from_start)/((double)window_size)));
}

/* Changed by Zhong Wang,3/3/2015
 *
 * Returns vector of counts for each windows in genomic data.
 *
 * To do this efficently, loop through the region once.
 *   Assume ... (1) user dosen't pass in something beyond the bounds of chrom_counts.
 *
 * Arguments:
 *  left_pos  --> the left position of each range request, absolute value, not relative to the read_from_bigWig_r window
 *  right_pos --> the right position of each range request, absolute value
 *  zoom      --> zoom_params_t
 *  chrom_counts --> raw_data_t, the results from read_from_bigWig_r
 *  dp        --> (return) genomic_data_point_t, output object.
 */
void get_genomic_data_2015(int left_pos, int right_pos, zoom_params_t zoom, raw_data_t chrom_counts, genomic_data_point_t dp) {
  init_genomic_data_point(dp, zoom);

  int left_idx  = left_pos - (chrom_counts.start + chrom_counts.offset);
  int right_idx = right_pos - (chrom_counts.start + chrom_counts.offset);

  // Sets up a bit of a  strange boundary condition, where 0's are included on all out-of-bounds windows.
  // After last night's experiment, I see this as preferable to throwing an error, however.
  if( right_idx >= chrom_counts.size) right_idx = chrom_counts.size;

  int center = floor((left_pos + right_pos)/2) - (chrom_counts.start + chrom_counts.offset);

  // Loop through incrementing each vector.
  for(int bp=left_idx; bp<right_idx; bp++) {

	// skip 0 reads. The performance is greatly improved.
	// Zhong Wang 6/20/2015
	// 1820 seconds  ==>170 seconds for 242078 * 441 (one CPU)
	if( (chrom_counts.forward[ bp ]==0 )  && (chrom_counts.reverse[ bp ]==0 ) )
		continue;

    for(int i=0;i<zoom.n_sizes;i++) {
      int which_bin= get_bin_number( center, bp, zoom.window_sizes[i], zoom.half_n_windows[i]);
      if(which_bin>=0 && bp>=0) {
        dp.forward[i][which_bin]+= (double)chrom_counts.forward[ bp ];
        dp.reverse[i][which_bin]+= (double)chrom_counts.reverse[ bp ];
      }
    }
  }

//Rprintf("==>get_genomic_data RD:(%d, %d), RANGE: (%d, %d), IDX (%d,%d),%d\n", chrom_counts.start, chrom_counts.offset, left_pos, right_pos, left_idx, right_idx, center);

}

// Improvement: Calculate range for each bin at each window scale.
// Zhong Wang 6/20/2015
// 1820 seconds  ==> 305  seconds for 242078 * 441 (one CPU)

void get_genomic_data_gpulike(int left_pos, int right_pos, zoom_params_t zoom, raw_data_t chrom_counts, genomic_data_point_t dp)
{
  init_genomic_data_point(dp, zoom);

  int center = floor((left_pos + right_pos)/2) - (chrom_counts.start + chrom_counts.offset);
  int right_idx = right_pos - (chrom_counts.start + chrom_counts.offset);
  if( right_idx >= chrom_counts.size) right_idx = chrom_counts.size;

  for(int i=0;i < zoom.n_sizes;i++) {
  for(int j=0;j < 2* zoom.half_n_windows[i];j++) {

	  int offset = 0;
	  if ( j>= zoom.half_n_windows[i] ) offset = 1;
	  int bin_start = center + zoom.window_sizes[i] * (-1*zoom.half_n_windows[i] + j) + offset;
	  int bin_stop  = center + zoom.window_sizes[i] * (-1*zoom.half_n_windows[i] + j + 1) + offset;

	  if( bin_start < 0)  bin_start = 0;
	  if( bin_stop >= right_idx )  bin_stop = right_idx;
	  if( bin_stop > bin_start)
  	  for(int k=bin_start; k<bin_stop; k++)
  	  {
        dp.forward[i][j] += (double)chrom_counts.forward[ k ];
        dp.reverse[i][j] += (double)chrom_counts.reverse[ k ];
      }
    }
  }
}

int** get_pre_bin( int max_dist, zoom_params_t zoom )
{
  int** pre_bin = (int**)Calloc(zoom.n_sizes, int* );

  for(int i=0;i<zoom.n_sizes;i++)
  {
	pre_bin[i] = (int*)Calloc( max_dist*2 + 1, int);
	memset( (void*)(pre_bin[i]), (int)-1, (size_t)(max_dist*2 + 1) );
  	for(int bp=0; bp<max_dist*2+1; bp++)
  	{
      pre_bin[i][bp] = get_bin_number( max_dist, bp, zoom.window_sizes[i], zoom.half_n_windows[i]);
  	}
  }

  return(pre_bin);
}

// Improvement: lookup the bin table pre-calculated for each position.
// Zhong Wang 6/20/2015
// 1820 seconds  ==>90 seconds for 242078 * 441 (one CPU)

void get_genomic_data_fast( int left_pos, int right_pos, zoom_params_t zoom, raw_data_t chrom_counts, genomic_data_point_t dp, int** pre_bin)
{
  init_genomic_data_point(dp, zoom);

  int left_idx  = left_pos - (chrom_counts.start ); // + chrom_counts.offset
  int right_idx = right_pos - (chrom_counts.start); //+ chrom_counts.offset

  // Sets up a bit of a  strange boundary condition, where 0's are included on all out-of-bounds windows.
  // After last night's experiment, I see this as preferable to throwing an error, however.
  if( right_idx >= chrom_counts.size) right_idx = chrom_counts.size;

  // Loop through incrementing each vector.
  for(int bp=left_idx;bp<right_idx; bp++) {
	if( bp<0 || ((chrom_counts.forward[ bp ]==0 )  && (chrom_counts.reverse[ bp ]==0 )) )
		continue;

    for(int i=0;i<zoom.n_sizes;i++) {
	  int which_bin = pre_bin[i][bp-left_idx];
      if(which_bin>=0) {
        dp.forward[i][which_bin]+= (double)chrom_counts.forward[ bp ];
        dp.reverse[i][which_bin]+= (double)chrom_counts.reverse[ bp ];
      }
    }
  }
}

/*
 * Returns the max value in data1 and data2.
 * WARNING: Use ONLY for positive integers 0 .. +Inf.
 */
double get_max(int n, double* data1, double* data2) {
  double max=-1;
  for (int i=0;i<n;i++) {
    if(data1[i] > max) max= data1[i];
	if(data2[i] > max) max= data2[i];
  }
  return(max);
}

double get_max_d1(int n, double* data1) {
  double max=-1;
  for (int i=0;i<n;i++) {
    if(data1[i] > max) max= data1[i];
  }
  if(max==0) max=1;
  else max= 0.05*max;
  return(max);
}

// CGD: 12-23-2013: Scale forward and reverse separately to avoid reverse signal at some highly paused genes?
void scale_genomic_data_strand_sep(zoom_params_t zoom, genomic_data_point_t dp) {
  double val_at_min=0.01;
  for(int i=0;i<zoom.n_sizes;i++) {
    // Get parameters.  Require value of 0.99 at MAX and 0.01 at 0.
    double max_val_fwd= get_max_d1(2*zoom.half_n_windows[i], dp.forward[i]);
    double max_val_rev= get_max_d1(2*zoom.half_n_windows[i], dp.reverse[i]);
    double alpha_fwd= /*2**/log(1/val_at_min - 1) / max_val_fwd; // Best without scaling here.
    double alpha_rev= /*2**/log(1/val_at_min - 1) / max_val_rev; // Best without scaling here.

	// Scale values with the logistic function.
    for(int j=0;j<2*zoom.half_n_windows[i]; j++) {
      dp.forward[i][j] = 1/ (1+ exp(-1*alpha_fwd*(dp.forward[i][j]-(max_val_fwd))));
      dp.reverse[i][j] = 1/ (1+ exp(-1*alpha_rev*(dp.reverse[i][j]-(max_val_rev))));
    }
  }
}


/*
 * Scales genomic data ... following a logistic function with the specified parameters...
 *  Pass this information in zoom_params_t(?!).  Or store it in genomic_data_point_t?  Or separate struct?!
 *
 * Scales using logistic function --> F(t)= 1/(1+e^(-\alpha(t-\beta)))
 *      Right now, \beta= MAX/2 (defines position of 0.5).
 *                 \alpha= log(1/0.01 - 1)/MAX  (Signal at 0 reads set to 0.01).
 *
 * 7/28/2013 -- Scaling best at 1*alpha and max set to 0.01 (AUC=0.9385945).
 * 		See 7/28/2013 notebook entry for all parameter settings tested.
 */
void scale_genomic_data_opt(zoom_params_t zoom, genomic_data_point_t dp) {
  double val_at_min=0.01;
  for(int i=0;i<zoom.n_sizes;i++) {
    // Get parameters.  Require value of 0.99 at MAX and 0.01 at 0.
    double max_val= get_max(2*zoom.half_n_windows[i], dp.forward[i], dp.reverse[i]);
	if(max_val == 0) max_val=1;
	else max_val= 0.05*max_val; // Try scaling to 5%.
    double alpha= /*2**/log(1/val_at_min - 1) / max_val; // Best without scaling here.

	// Scale values with the logistic function.
    for(int j=0;j<2*zoom.half_n_windows[i]; j++) {
      dp.forward[i][j] = 1/ (1+ exp(-1*alpha*(dp.forward[i][j]-(max_val))));
      dp.reverse[i][j] = 1/ (1+ exp(-1*alpha*(dp.reverse[i][j]-(max_val))));
    }
  }
}


/*
  Note: Comparing SVM optimized and true logistic.  The following R shows off the difference between the two
  approaches:

	scaled_logistic_function <- function(x) {
	  MAX <- max(x)
	  beta_ <- MAX*0.5
	  alpha_ <- log(1/0.01 - 1)/beta_
	  1/(1+ exp(-1*alpha_*(x-beta_))) ## 2 is an arbitrarily chosen value...
	}
	scaled_logistic_function(c(1:10))

	scaled_logistic_function_ <- function(x) {
	  MAX <- max(x)
	  beta_ <- MAX*0.5
	  alpha_ <- log(1/0.01 - 1)/beta_
	  1/(1+ exp(-1*alpha_*(x-beta_/2))) ## 2 is an arbitrarily chosen value...
	}
	scaled_logistic_function_(c(1:10))-scaled_logistic_function(c(1:10))
	plot(scaled_logistic_function(c(1:10)), type="b")
	points(scaled_logistic_function_(c(1:10)), type="b", col="dark red")

  It's a complete left shift, but a bit straing in that the bottom value is no longer 0.  Try shifting further left using MAX*0.2?!
 */

void scale_genomic_data_simple_max(zoom_params_t zoom, genomic_data_point_t dp) {
  for(int i=0;i<zoom.n_sizes;i++) {
    double max_val= get_max(2*zoom.half_n_windows[i], dp.forward[i], dp.reverse[i]);
    if(max_val == 0) max_val=1;

    for(int j=0;j<2*zoom.half_n_windows[i]; j++) {
      dp.forward[i][j] = dp.forward[i][j]/max_val;
      dp.reverse[i][j] = dp.reverse[i][j]/max_val;
    }
  }
}

/*
 * Moves C genomic_data_point_t type to a SEXP for return to R.
 */
SEXP data_point_to_r_list(zoom_params_t zoom, genomic_data_point_t dp) {
  SEXP data_point;
  protect(data_point = allocVector(VECSXP, 2*zoom.n_sizes));

  for(int i=0;i<zoom.n_sizes;i++) {
    // Create R object.
    SEXP size_t_for, size_t_rev;
    protect(size_t_for = allocVector(REALSXP, zoom.half_n_windows[i]*2));
    protect(size_t_rev = allocVector(REALSXP, zoom.half_n_windows[i]*2));
    SET_VECTOR_ELT(data_point, 2*i, size_t_for);
    SET_VECTOR_ELT(data_point, 2*i+1, size_t_rev);

    // Copy data from dp to R object.
    double *size_t_for_c = REAL(size_t_for);
    double *size_t_rev_c = REAL(size_t_rev);
    for(int j=0;j<2*zoom.half_n_windows[i];j++) {
      size_t_for_c[j] = dp.forward[i][j];
      size_t_rev_c[j] = dp.reverse[i][j];
    }
	UNPROTECT(2);
  }
  UNPROTECT(1);
  return(data_point);
}

/*
 * Moves C genomic_data_point_t type to a SEXP for return to R.
 * This function generates a single vector.
 */
SEXP data_point_to_r_vect(zoom_params_t zoom, genomic_data_point_t dp) {
  SEXP data_point;

  // Count  number of windows to allocate R vector...
  int n_windows=0;
  for(int i=0;i<zoom.n_sizes;i++)
    n_windows+= zoom.half_n_windows[i]*2; // *2 for full set of windows.

  protect(data_point = allocVector(REALSXP, 2*n_windows)); // *2 for both strands.
  double *data_point_c = REAL(data_point);

  int k=0;
  for(int i=0;i<zoom.n_sizes;i++) {
    for(int j=0;j<2*zoom.half_n_windows[i];j++) {
      data_point_c[k] = dp.forward[i][j];
      data_point_c[n_windows+k++] = dp.reverse[i][j];
    }
  }
  return(data_point);
}

/*
 * Reads the specified region from a bigWig file.
 */
raw_data_t read_from_bigWig_r(const char *chrom, int start, int end, bigWig_t *bw_fwd, bigWig_t *bw_rev) {
  raw_data_t rd;
  int out_length, out_is_blank;

  rd.start = start;
  if(start < 0) {
    rd.offset= -1*(start);
	rd.start=0;
  } else {
    rd.offset= 0;
  }

  // Read in raw data from the bigWig.
  rd.forward= bigwig_readi(bw_fwd, chrom, start, end, 1, 1, &out_length, &out_is_blank);
  rd.reverse= bigwig_readi(bw_rev, chrom, start, end, 1, 1, &out_length, &out_is_blank);
  rd.size= out_length;

  return(rd);
}

/* Bigwiglib allocs from the C memory pool. */
void free_raw_data(raw_data_t rd) {
  free(rd.forward);
  free(rd.reverse);
}

/* Added by Zhong Wang 3/3/2015

	Merge multiple ranges into one big range to improve the efficiency of bigWig reading

	SEXP chrom_r      --> the chromosome part of allranges transfered by R
	SEXP centers_r    --> the center position part of allranges transfered by R
	int max_dist      --> range width
	const char* pszChr--> current chromosome no.
	int nrange_center --> current range center
	int nrange_max    --> the biggest width for one time bigWig reading
	int nLocal_id     --> current group ID to identify group number of each requested range and indicate how many times have been called to this function.
	int* pLocal_start --> (return)suggested new start position of  merged big range
	int* pLocal_stop  --> (return)suggested new stop position of  merged big range
	int* pLocal_range_table--> (return)the group table for each range, 0 means no group associated with this range.
*/
int merge_adjacent_range( SEXP chrom_r, SEXP centers_r, int* pStart_idx, int max_dist, const char* pszChr, int nrange_center, int nrange_max, int nLocal_id, int* pLocal_start, int* pLocal_stop, int* pLocal_range_table )
{
    int n_centers = Rf_nrows(centers_r);
    int *centers = INTEGER(centers_r);

	int start_pos = nrange_center;
	int stop_pos = nrange_center;
	int nRange = 0;
	for(int i = *pStart_idx; i<n_centers; i++ )
	{
		// The range has been retrived, skip it.
		if( pLocal_range_table[i] !=0 ) continue;

    	const char* pszChr_i = CHAR(STRING_ELT(chrom_r,i));
		// not at same chromosom, skip it.
		if( strcmp(pszChr_i, pszChr) != 0) break;

		// if within the max range at one time reading
		if( nrange_center - nrange_max <= centers[i] - max_dist  && centers[i] + max_dist <= nrange_center + nrange_max)
		{
			pLocal_range_table[i] = nLocal_id;

			if( centers[i] - max_dist < start_pos ) start_pos = centers[i] - max_dist;
			if( centers[i] + max_dist > stop_pos )  stop_pos  = centers[i] + max_dist;
			nRange++;
			*pStart_idx = i;

			if(nRange>5000) break;
		}
		else
			break;
	}

	*pLocal_start = start_pos;
	*pLocal_stop = stop_pos;

// Rprintf("==>LAST %d =(%d,%d)\n", nLocal_id, *pLocal_start, *pLocal_stop);

	return(nRange);
}

/* Changed by Zhong Wang 3/3/2015
 *
 * R entry point ... for getting a particular center (or vector of centers).
 *
 * Switch to R vector using:
 * t(matrix(unlist(list(c(1:10), c(11:20), c(0:9))), ncol=3))
 */
SEXP get_genomic_data_R(SEXP chrom_r, SEXP centers_r, SEXP bigwig_plus_file_r, SEXP bigwig_minus_file_r, SEXP model_r, SEXP scale_r) {
  bool bscale= BOOLEAN_ELT(scale_r, 0);
  int n_centers = Rf_nrows(centers_r);
  int *centers = INTEGER(centers_r);

  // Set up model variable.
  zoom_params_t zoom;
  zoom.n_sizes= Rf_nrows(VECTOR_ELT(model_r, 0));
  zoom.window_sizes= INTEGER(VECTOR_ELT(model_r, 0));
  zoom.half_n_windows= INTEGER(VECTOR_ELT(model_r, 1));

  // Open bigWig files.
  //
  assert(is_bigwig(CHAR(STRING_ELT(bigwig_plus_file_r, 0)))==1 && is_bigwig(CHAR(STRING_ELT(bigwig_minus_file_r, 0)))==1);

  bigWig_t *bw_fwd = bigwig_load(CHAR(STRING_ELT(bigwig_plus_file_r,  0)), ".");
  bigWig_t *bw_rev = bigwig_load(CHAR(STRING_ELT(bigwig_minus_file_r, 0)), ".");

  if (bw_fwd==NULL || bw_rev==NULL) return( R_NilValue );

  int max_dist= max_dist_from_center(zoom.n_sizes, zoom.window_sizes, zoom.half_n_windows);
  int** pre_bin = get_pre_bin( max_dist, zoom );

  // Set up return variable.
  genomic_data_point_t dp= alloc_genomic_data_point(zoom);

  //SEXP processed_data = PROTECT(allocVector(VECSXP, n_centers));
  int n_windows=0;
  for(int i=0;i<zoom.n_sizes;i++)
    n_windows+= zoom.half_n_windows[i]*2; // *2 (half*2) * 2(minus, plus)
  SEXP processed_data = PROTECT(allocMatrix(REALSXP, n_windows*2, n_centers ));
  double *pRMat = REAL(processed_data);

  int* pLocal_range_table = Calloc(n_centers, int);

  int nLocal_id=1;
  int k=0;
  while( k < n_centers )
  {
	// group_id>0, this range has been merged.
	if( pLocal_range_table[k]>0 ) continue;

    const char* pszChr = CHAR(STRING_ELT(chrom_r,k));
    int nLocal_start =0;
    int nLocal_stop =0;
    // The biggest range allowed by programmer; 1M, maybe too small?
	int nMerge_max = MERGE_RANGE;
	if( nMerge_max < max_dist*8 ) nMerge_max = max_dist*8;

	int start_idx = k;
	int nRangeCnt = merge_adjacent_range( chrom_r, centers_r, &start_idx, max_dist, pszChr, centers[k], nMerge_max, nLocal_id, &nLocal_start, &nLocal_stop, pLocal_range_table);
    if( nRangeCnt<=0)
    {
		// Although impossible, but we need to watch out.
		Rprintf("\nAn error occurs in the function of get_genomic_data_R, please contact the author!\n");
		Free(pLocal_range_table);
		return(R_NilValue);
	}

	if(nLocal_start<0) nLocal_start=0;

	// Read bigwig using a big merged range
    raw_data_t rd_local= read_from_bigWig_r( pszChr, nLocal_start, nLocal_stop, bw_fwd, bw_rev);

    for(int i=k;i<=start_idx;i++)
    {
	  // if not the curent id, done or no group, skip the range
	  if(pLocal_range_table[i] != nLocal_id) continue;

	  // Read raw data, do windowing specified in model_r, and scale.

	  //get_genomic_data_2015( centers[i] - max_dist, centers[i] + max_dist, zoom, rd_local, dp ); // Data center should be @ max_dist.  Total window should be 2*max_dist.  Center relative to read_from_bigWig_r, rather than chromosome.
	  //get_genomic_data_gpulike( centers[i] - max_dist, centers[i] + max_dist, zoom, rd_local, dp ); // Data center should be @ max_dist.  Total window should be 2*max_dist.  Center relative to read_from_bigWig_r, rather than chromosome.
      get_genomic_data_fast( centers[i] - max_dist, centers[i] + max_dist, zoom, rd_local, dp, pre_bin );

      if(bscale) scale_genomic_data_strand_sep( zoom, dp); //scale_genomic_data_opt

	  // Record ...
	  //SEXP data_point= data_point_to_r_vect( zoom, dp );
	  //SET_VECTOR_ELT(processed_data, i, data_point);
	  //UNPROTECT(1);
      int zk=0;
      int base = i*n_windows*2;
      for(int zi=0; zi<zoom.n_sizes; zi++)
      for(int zj=0; zj<2*zoom.half_n_windows[zi]; zj++) {
        pRMat[base + zk] = dp.forward[zi][zj];
        pRMat[base + n_windows + zk++] = dp.reverse[zi][zj];
      }
    }

    free_raw_data(rd_local);

    nLocal_id++;
    k = start_idx + 1;
  }

  //Calloc and Free are paired functions
  Free(pLocal_range_table);

  for(int i=0;i<zoom.n_sizes;i++) Free(pre_bin[i]);
  Free(pre_bin);

  free_genomic_data_point(dp, zoom);
  bigwig_free(bw_fwd);
  bigwig_free(bw_rev);
  UNPROTECT(1);

  return(processed_data);
}

/*
 * We need a simple vector to pass into the more general (??) functions.
 */
/*int *genomic_data_point_to_vector() {

  int *c_list = (int*)calloc((2*zoom.n_prev_bins[zoom.n_sizes]), sizeof(int)); // Because I'm going to destroy it ...

}*/
