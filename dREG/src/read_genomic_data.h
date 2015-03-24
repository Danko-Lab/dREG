#ifndef read_genomic_data_H
#define read_genomic_data_H
#include "bigwiglib.h"

typedef struct {
  int n_sizes;
  int *window_sizes;
  int *half_n_windows;

  int *n_prev_bins;
} zoom_params_t;

typedef struct {
  double **forward;
  double **reverse;
} genomic_data_point_t;

typedef struct {
  int start;
  int size;
  int offset;  // If raw_data_t cannot start at the left endge of the window, this value indicates the start position.

  int *forward;
  int *reverse;
} raw_data_t;

raw_data_t read_from_bigWig_r(const char *chrom, int start, int end, bigWig_t *bw_fwd, bigWig_t *bw_rev);
void get_genomic_data(int left_pos, int right_pos, zoom_params_t zoom, raw_data_t chrom_counts, genomic_data_point_t dp);
void scale_genomic_data(zoom_params_t zoom, genomic_data_point_t dp);
int max_dist_from_center(int n_sizes, int *window_sizes, int *half_n_windows);
SEXP get_genomic_data_R(SEXP chrom_r, SEXP centers_r, SEXP bigwig_plus_file_r, SEXP bigwig_minus_file_r, SEXP model_r);

#endif

