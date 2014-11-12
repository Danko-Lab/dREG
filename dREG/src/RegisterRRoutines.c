/**************************************************************
 *
 *      Associates a vector of genomic featuers (e.g. genes, CpG islands, etc.) with a table of sequence reads.
 *
 **************************************************************/
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "read_genomic_data.h"

/**************************************************************
 *
 *      Register entry points...
 *
 **************************************************************/
void R_init_featureDetector(DllInfo *info) {
     R_CallMethodDef callMethods[]  = {
       {"get_genomic_data_R", (DL_FUNC)&get_genomic_data_R, 5},
       {NULL, NULL, 0}
     };

        R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
