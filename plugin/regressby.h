#ifndef REGRESSBY_H
#define REGRESSBY_H    1

#include <locale.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>

#include "regresstools.c"

// Container structure for Stata-provided info
struct StataInfo {
    GT_size in1;
    GT_size in2;
    GT_size N;
    GT_size J;
    GT_size nj_max;
    GT_size free;
    //
    GT_bool constant;
    GT_bool covariance;
    GT_bool robust;
    GT_bool cluster;
    GT_bool benchmark;
    GT_size kx;
    //
    //
    GT_size *nobs;
    GT_size *info;
    GT_size *cinfo;
    GT_size *groupvar;
    GT_size *clustervar;
    //
    gsl_matrix *X;
    gsl_vector *y;
};

// Main functions
ST_retcode ssf_parse_info    (struct StataInfo *st_info);
ST_retcode ssf_read_varlist  (struct StataInfo *st_info);
ST_retcode ssf_regressby     (struct StataInfo *st_info, char *fname);
void ssf_free (struct StataInfo *st_info);

#endif /* regressby.h */
