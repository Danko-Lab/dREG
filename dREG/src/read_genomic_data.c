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
 * get_bin_number --> Returns the bin that a particular posotion falls into.
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

/*
 * Returns vector of counts for each windows in genomic data.
 *
 * To do this efficently, loop through the region once.  
 *   Assume ... (1) user dosen't pass in something beyond the bounds of chrom_counts.
 *
 * Arguments: 
 *  center  --> Center is relative to the read_from_bigWig_r window, rather than the chromosome.
 *  n_sizes --> Number of 
 */
void get_genomic_data(int center, zoom_params_t zoom, raw_data_t chrom_counts, genomic_data_point_t dp) {
  init_genomic_data_point(dp, zoom);

  // Get the max boundary of our window.
  int max_bounds = max_dist_from_center(zoom.n_sizes, zoom.window_sizes, zoom.half_n_windows);
  int left_edge= center - max_bounds + chrom_counts.offset;
  int right_edge= center + max_bounds+1;
  
  // Sets up a bit of a  strange boundary condition, where 0's are included on all out-of-bounds windows.
  // After last night's experiment, I see this as preferable to throwing an error, however.
  if(right_edge >= chrom_counts.size) right_edge = chrom_counts.size;

  // Loop through incrementing each vector.
  for(int bp= left_edge;bp<right_edge;bp++) {
    for(int i=0;i<zoom.n_sizes;i++) {
      int which_bin= get_bin_number(center, bp, zoom.window_sizes[i], zoom.half_n_windows[i]);
      if(which_bin>=0) {
        dp.forward[i][which_bin]+= (double)chrom_counts.forward[bp-chrom_counts.offset];
        dp.reverse[i][which_bin]+= (double)chrom_counts.reverse[bp-chrom_counts.offset];
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

  It's a complete left shift, but a bit straing in that the bottom value is no loger 0.  Try shifting further left using MAX*0.2?!
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
    // Creat R object.
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
  
  if(start < 0) {
    rd.offset= -1*(start);
	start=0;
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

/*
 * R entry point ... for getting a particular center (or vector of centers).
 *
 * Switch to R vector using:
 * t(matrix(unlist(list(c(1:10), c(11:20), c(0:9))), ncol=3))
 */
SEXP get_genomic_data_R(SEXP chrom_r, SEXP centers_r, SEXP bigwig_plus_file_r, SEXP bigwig_minus_file_r, SEXP model_r) {
  int n_centers = Rf_nrows(centers_r);
  int *centers = INTEGER(centers_r);
  
  // Set up model variable.
  zoom_params_t zoom;
  zoom.n_sizes= Rf_nrows(VECTOR_ELT(model_r, 0));
  zoom.window_sizes= INTEGER(VECTOR_ELT(model_r, 0));
  zoom.half_n_windows= INTEGER(VECTOR_ELT(model_r, 1));

  // Open bigWig files.
  assert(is_bigwig(CHAR(STRING_ELT(bigwig_plus_file_r, 0)))==1 && is_bigwig(CHAR(STRING_ELT(bigwig_minus_file_r, 0)))==1);
  bigWig_t *bw_fwd = bigwig_load(CHAR(STRING_ELT(bigwig_plus_file_r,  0)), ".");
  bigWig_t *bw_rev = bigwig_load(CHAR(STRING_ELT(bigwig_minus_file_r, 0)), ".");
  
  // Set up return variable.
  genomic_data_point_t dp= alloc_genomic_data_point(zoom);
  SEXP processed_data;
  PROTECT(processed_data = allocVector(VECSXP, n_centers));

  for(int i=0;i<n_centers;i++) {
    // Read raw data, do windowing specified in model_r, and scale.
    int max_dist= max_dist_from_center(zoom.n_sizes, zoom.window_sizes, zoom.half_n_windows);
    raw_data_t rd= read_from_bigWig_r(CHAR(STRING_ELT(chrom_r,i)), centers[i]-max_dist, centers[i]+max_dist, bw_fwd, bw_rev);
    get_genomic_data(max_dist, zoom, rd, dp); // Data center should be @ max_dist.  Total window should be 2*max_dist.  Center relative to read_from_bigWig_r, rather than chromosome.
    free_raw_data(rd);
    scale_genomic_data_strand_sep(zoom, dp); //scale_genomic_data_opt

    // Record ...
    SEXP data_point= data_point_to_r_vect(zoom, dp);//data_point_to_list(zoom, dp);
    SET_VECTOR_ELT(processed_data, i, data_point);
    UNPROTECT(1);
  }
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
