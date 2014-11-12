#include "bigwiglib.h"

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "udc.h"
#include "bigWig.h"
#include "obscure.h"
#include "errCatch.h"
#include "hmmstats.h"
#include "localmem.h"

int is_bigwig(const char * filename) {
  return isBigWig((char*)filename);
}

bigWig_t * bigwig_load(const char * filename, const char * udc_dir) {
  bigWig_t * bigwig;
  struct errCatch * err;

  /* set cache */
  if (udc_dir != NULL)
    udcSetDefaultDir((char*) udc_dir);

  /* setup error management & try to open file */
  err = errCatchNew();
  if (errCatchStart(err))
    bigwig = bigWigFileOpen((char*)filename);
  errCatchEnd(err);
  if (err->gotError) {
    fprintf(stderr, "error: %s\n", err->message->string);
    errCatchFree(&err);
    return NULL;
  }
  errCatchFree(&err);

  return bigwig;
}

void bigwig_free(bigWig_t * bw) {
  if (bw != NULL)
    bbiFileClose(&bw);
}

int bigwig_valid_chrom(bigWig_t * bw, char * chrom) {
  struct bbiChromInfo * chrom_i, * chromList = bbiChromList(bw);
  int result = 0;

  for (chrom_i = chromList; chrom_i != NULL; chrom_i = chrom_i->next)
    if (!strcmp(chrom_i->name, chrom)) {
      result = 1;
      break;
    }

  bbiChromInfoFreeList(&chromList);

  return result;
}

double * bigwig_readf(bigWig_t * bw, char * chrom, int start, int end, int step, int abs, int * out_length, int * out_is_blank) {
  struct lm * localMem = lmInit(0); /* use default value */
  struct bbiInterval * intervals 
    = bigWigIntervalQuery(bw, chrom, start, end, localMem);
  int nIntervals = slCount(intervals);
  int size = (end - start)/step;
  double * result = (double*) calloc(size, sizeof(double));

  if (nIntervals > 0) {
    struct bbiInterval * interval;
    int count = 0;
    double sum = 0.0;
    double left, right;
    int idx;
    
    left = start;
    right = start + step;
    idx = 0;

    for (interval = intervals; interval != NULL && idx < size; interval = interval->next) {
      /* interval starts beyond current step */
      if (((double)interval->start) >= right) {
	
	/* save current value */
	if (count > 0)
	  result[idx] = sum;
	
	count = 0;
	sum = 0.0;
	while (idx < size && ((double)interval->start) >= right) {
	  ++idx;
	  left += step;
	  right += step;
	}
      }

      /* interval starts at or before the current step */
      if (((double)interval->start) < right) {
	sum += interval->val;
	++count;
      }

      /* interval ends beyond the current step */
      if (((double)interval->end) > right && idx < size) {
	do {
	  /* save current step */
	  if (count > 0)
	    result[idx] = sum;
	    
	  ++idx;
	  left += step;
	  right += step;
	    
	  count = (((double) interval->end > right) ? 1 : 0);
	  sum = interval->val;
	} while (idx < size && ((double)interval->end > right));
      }
    }
    
    if (count > 0 && idx < size)
      result[idx] = sum;

    /* do abs */
    if (abs == 1) {
      int i;
      
      for (i = 0; i < size; ++i)
	result[i] = fabs(result[i]);
    }

    *out_is_blank = 0;
  } else
    *out_is_blank = 1;


  /* clean-up */
  lmCleanup(&localMem);

  *out_length = size;
  return result;
}

int * bigwig_readi(bigWig_t * bw, char * chrom, int start, int end, int step, int abs, int * out_length, int * out_is_blank) {
  double * data = bigwig_readf(bw, chrom, start, end, step, abs, out_length, out_is_blank);
  int * result;
  int i, len;

  if (data == NULL)
    return NULL;
  
  len = *out_length;
  result = (int*) calloc(len, sizeof(int));

  if (*out_is_blank == 0) {
    for (i = 0; i < len; ++i)
      result[i] = (int) round(data[i]);
  }

  /* clean-up */
  free(data);

  return result;
}
