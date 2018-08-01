/*********************************************************************
 * Program: regressby.c
 * Blame:   Mauricio Caceres Bravo <mauricio.caceres.bravo@gmail.com>
 * Created: Tue Jul 31 21:11:03 EDT 2018
 * Updated: Wed Aug  1 09:12:03 EDT 2018
 * Purpose: Stata plugin for fast regression by group
 * Note:    See stata.com/plugins for more on Stata plugins
 * Version: 0.1.0
 *********************************************************************/

/**
 * @file regressby.c
 * @author Mauricio Caceres Bravo
 * @date 1 Aug 2018
 * @brief Stata plugin
 *
 * This file should only ever be called from regressby.ado
 *
 * @see help regressby
 * @see http://www.stata.com/plugins for more on Stata plugins
 */

#include "regressby.h"

int main()
{
    return(0);
}

int WinMain()
{
    return(0);
}

STDLL stata_call(int argc, char *argv[])
{
    ST_retcode rc = 0;
    setlocale(LC_ALL, "");
    struct StataInfo *st_info = malloc(sizeof(*st_info));
    st_info->free = 0;

    clock_t timer  = clock();
    size_t flength = strlen(argv[0]) + 1;
    char *fname    = malloc(sizeof(char) * flength);
    memset (fname, '\0', sizeof(char) * flength);
    strcpy (fname, argv[0]);

    if ( (rc = ssf_parse_info (st_info)) ) {
        goto exit;
    }

    if ( st_info->benchmark )
        sf_running_timer (&timer, "\tPlugin step 1: Parse stata info");

    if ( (rc = ssf_read_varlist (st_info)) ) {
        goto exit;
    }

    if ( st_info->benchmark )
        sf_running_timer (&timer, "\tPlugin step 2: Read in variables");

    if ( (rc = ssf_regressby (st_info, fname)) ) {
        goto exit;
    }

    if ( st_info->benchmark )
        sf_running_timer (&timer, "\tPlugin step 3: Ran regressions");

exit:
    ssf_free (st_info);
    free (st_info);
    return (rc);
}

/**
 * @brief Parse variable info from Stata
 *
 * @param st_info Pointer to container structure for Stata info
 * @return Stores in @st_info various info from Stata for the pugin run
 */
ST_retcode ssf_parse_info (struct StataInfo *st_info)
{
    ST_retcode rc = 0;
    GT_size in1, in2, N;
    GT_size constant,
            covariance,
            robust,
            cluster,
            benchmark,
            kx;

    in1 = SF_in1();
    in2 = SF_in2();
    N   = in2 - in1 + 1;

    if ( (rc = sf_scalar_size("regressby_constant",   &constant)   )) goto exit;
    if ( (rc = sf_scalar_size("regressby_covariance", &covariance) )) goto exit;
    if ( (rc = sf_scalar_size("regressby_robust",     &robust)     )) goto exit;
    if ( (rc = sf_scalar_size("regressby_cluster",    &cluster)    )) goto exit;
    if ( (rc = sf_scalar_size("regressby_benchmark",  &benchmark)  )) goto exit;
    if ( (rc = sf_scalar_size("regressby_kx",         &kx)         )) goto exit;

    st_info->in1 = in1;
    st_info->in2 = in2;
    st_info->N   = N;

    st_info->constant   = constant;
    st_info->covariance = covariance;
    st_info->robust     = robust;
    st_info->cluster    = cluster;
    st_info->benchmark  = benchmark;
    st_info->kx         = kx;

exit:
    return(rc);
}

/**
 * @brief Read varlist to hash from stata
 *
 * @param st_info Pointer to container structure for Stata info
 * @return Stores in @st_info->st_charx the input variables
 */
ST_retcode ssf_read_varlist (struct StataInfo *st_info)
{
    ST_retcode rc = 0;
    ST_double z;
    GT_size i, j, k, sel, *nptr;

    // Allocate memory for regression variables
    st_info->groupvar   = calloc(st_info->N, sizeof st_info->groupvar);
    st_info->clustervar = calloc(st_info->cluster? st_info->N: 1, sizeof st_info->clustervar);

    st_info->X = gsl_matrix_alloc (st_info->N, st_info->kx);
    st_info->y = gsl_vector_alloc (st_info->N);

    if ( st_info->groupvar   == NULL ) return (sf_oom_error("ssf_read_varlist", "groupvar"));
    if ( st_info->clustervar == NULL ) return (sf_oom_error("ssf_read_varlist", "clustervar"));

    if ( st_info->X == NULL ) return (sf_oom_error("ssf_read_varlist", "X"));
    if ( st_info->y == NULL ) return (sf_oom_error("ssf_read_varlist", "y"));

    st_info->free = 1;

    // Read in variables
    if ( st_info->cluster ) {
        for (i = 0; i < st_info->N; i++) {
            sel = i + st_info->in1;

            // First variable is group
            if ( (rc = SF_vdata(1, sel, &z)) ) goto exit;
            st_info->groupvar[i] = (GT_size) z;

            // Second variable is cluster
            if ( (rc = SF_vdata(2, sel, &z)) ) goto exit;
            st_info->clustervar[i] = (GT_size) z;

            // Third variable is y
            if ( (rc = SF_vdata(3, sel, &z)) ) goto exit;
            gsl_vector_set (st_info->y, i, z);

            // Last kx variables are covariates
            for (k = 0; k < st_info->kx; k++) {
                if ( (rc = SF_vdata(4 + k, sel, &z)) ) goto exit;
                gsl_matrix_set (st_info->X, i, k, z);
            }
        }
    }
    else {
        for (i = 0; i < st_info->N; i++) {
            sel = i + st_info->in1;

            // First variable is group
            if ( (rc = SF_vdata(1, sel, &z)) ) goto exit;
            st_info->groupvar[i] = (GT_size) z;

            // Second variable is y
            if ( (rc = SF_vdata(2, sel, &z)) ) goto exit;
            gsl_vector_set (st_info->y, i, z);

            // Last kx variables are covariates
            for (k = 0; k < st_info->kx; k++) {
                if ( (rc = SF_vdata(3 + k, sel, &z)) ) goto exit;
                gsl_matrix_set (st_info->X, i, k, z);
            }
        }
    }

    // Allocate memory for panel info; note that the last observation of
    // groupvar is J by construction
    st_info->J    = st_info->groupvar[st_info->N - 1];
    st_info->info = calloc(st_info->J + 1, sizeof st_info->info);
    st_info->nobs = calloc(st_info->J,     sizeof st_info->nobs);

    if ( st_info->info == NULL ) return (sf_oom_error("ssf_read_varlist", "info"));
    if ( st_info->nobs == NULL ) return (sf_oom_error("ssf_read_varlist", "nobs"));

    st_info->free = 2;

    // Panel setup from groupvar
    st_info->J = sf_panelsetup (st_info->groupvar, st_info->info, st_info->N);

    // Largest group
    nptr = st_info->nobs;
    st_info->nj_max = 0;
    for (j = 0; j < st_info->J; j++, nptr++) {
        *nptr = st_info->info[j + 1] - st_info->info[j];
        if ( *nptr > st_info->nj_max ) st_info->nj_max = *nptr;
    }

    // Space for cluster info, should there be any
    st_info->cinfo = calloc(st_info->cluster? st_info->nj_max + 1: 1, sizeof st_info->cinfo);
    if ( st_info->cinfo == NULL ) return (sf_oom_error("ssf_read_varlist", "cinfo"));
    st_info->free = 3;

exit:
    return(rc);
}

/**
 * @brief Read varlist to hash from stata
 *
 * @param st_info Pointer to container structure for Stata info
 * @return Stores in @st_info->st_charx the input variables
 */
ST_retcode ssf_regressby (struct StataInfo *st_info, char *fname)
{
    ST_retcode rc = 0;
    ST_double z, *cptr, J_dbl, Jout_dbl, kout_dbl;
    GT_size i, j, k, nj, nc, kout, Jout = 0;

    GT_bool hac0 = (st_info->robust == 0) & (st_info->cluster == 0);
    GT_bool hac1 = (st_info->robust == 1) & (st_info->cluster == 0);

    gsl_matrix_view X;
    gsl_vector_view y;
    gsl_matrix_view Xi;

    gsl_vector *Xy = gsl_vector_alloc (st_info->kx);
    gsl_vector *b  = gsl_vector_alloc (st_info->kx);

    gsl_vector_set_zero (b);
    gsl_vector_set_zero (Xy);

    gsl_matrix *Xe    = gsl_matrix_alloc (1,           st_info->kx);
    gsl_matrix *A     = gsl_matrix_alloc (st_info->kx, st_info->kx);
    gsl_matrix *XeX   = gsl_matrix_alloc (hac1? st_info->kx: 1, hac1? st_info->kx: 1);
    gsl_matrix *Vaux  = gsl_matrix_alloc (hac1? st_info->kx: 1, hac1? st_info->kx: 1);
    gsl_matrix *Vb;

    gsl_matrix_set_zero (A);
    gsl_matrix_set_zero (XeX);
    gsl_matrix_set_zero (Vaux);

    kout  = 2 * st_info->kx;
    kout += st_info->covariance? (st_info->kx * (st_info->kx - 1) / 2): 0;
    ST_double *coefs = calloc (st_info->J * kout, sizeof(coefs));

    cptr = coefs;
    for (j = 0; j < st_info->J; j++) {
        nj = st_info->nobs[j];
        if ( nj < st_info->kx )
            continue;

        // Views of y, X for the jth group
        X = gsl_matrix_submatrix(
            st_info->X,
            st_info->info[j], 0,
            nj, st_info->kx
        );

        y = gsl_vector_subvector(
            st_info->y,
            st_info->info[j],
            nj
        );

        // Set A = X' X and Xy = X' y
        if ( (rc = gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &X.matrix, &X.matrix, 0.0, A)) ) {
            goto exit;
        }
        if ( (rc = gsl_blas_dgemv(CblasTrans, 1.0, &X.matrix, &y.vector, 0.0, Xy)) ) {
            goto exit;
        }

        // Cholesky decomposition; b will contain the OLS coefficients
        // and A will contain (X' X)^-1

        if ( (rc = gsl_linalg_cholesky_decomp1(A))      ) goto exit;
        if ( (rc = gsl_linalg_cholesky_solve(A, Xy, b)) ) goto exit;
        if ( (rc = gsl_linalg_cholesky_invert(A))       ) goto exit;

        // Save the coefficients
        for (k = 0; k < st_info->kx; k++, cptr++) {
            *cptr = gsl_vector_get(b, k);
        }

        // y contains the errors
        if ( (rc = gsl_blas_dgemv (CblasNoTrans, -1.0, &X.matrix, b, 1.0, &y.vector)) ) {
            goto exit;
        }

        // Compute vcov matrix; Vb will point to V(b)
        if ( hac0 ) {
            // sigma is (e' e) / (n - k)
            if ( (rc = gsl_blas_ddot (&y.vector, &y.vector, &z)) ) goto exit;
            if ( (rc = gsl_matrix_scale (A, z / ((ST_double) nj - st_info->kx))) ) goto exit;
            Vb = A;
        }
        else if ( hac1 ) {
            // Set X to  e :* X (each row times e)
            for (i = 0; i < nj; i++) {
                z  = gsl_vector_get(&y.vector, i);
                Xi = gsl_matrix_submatrix(
                    &X.matrix,
                    i, 0,
                    1, st_info->kx
                );
                if ( (rc = gsl_matrix_scale (&Xi.matrix, z)) ) goto exit;
            }
            if ( (rc = gsl_blas_dgemm (CblasTrans,   CblasNoTrans, 1.0, &X.matrix, &X.matrix, 0.0, XeX))  ) goto exit;
            if ( (rc = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, A,         XeX,       0.0, Vaux)) ) goto exit;
            if ( (rc = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Vaux,      A,         0.0, XeX))  ) goto exit;
            if ( (rc = gsl_matrix_scale (XeX, (nj / ((ST_double) nj - st_info->kx)))) ) goto exit;
            Vb = XeX;
        }
        else {
            nc = sf_panelsetup(st_info->clustervar + st_info->info[j], st_info->cinfo, nj);
            nc++;
            nc--;
            Vb = A;
            // TODO: Cluster SE
        }

        // Save se
        for (k = 0; k < st_info->kx; k++, cptr++) {
            *cptr = sqrt(gsl_matrix_get(Vb, k, k));
        }

        // Save covariance
        for (i = 0; i < st_info->kx; i++) {
            for (k = 0; k < i; k++, cptr++) {
                *cptr = gsl_matrix_get(Vb, i, k);
            }
        }

        Jout++;
    }

    // Save output
    FILE *fhandle = fopen(fname, "wb");

    J_dbl    = (ST_double) st_info->J;
    Jout_dbl = (ST_double) Jout;
    kout_dbl = (ST_double) kout;

    fwrite (&J_dbl, sizeof(J_dbl), 1, fhandle);
    for (j = 0; j < st_info->J; j++) {
        z = (ST_double) st_info->nobs[j];
        fwrite (&z, sizeof(z), 1, fhandle);
    }

    fwrite (&Jout_dbl, sizeof(Jout_dbl), 1, fhandle);
    fwrite (&kout_dbl, sizeof(kout_dbl), 1, fhandle);
    fwrite (coefs, sizeof(coefs), Jout * kout, fhandle);

    fclose (fhandle);

exit:
    gsl_vector_free (Xy);
    gsl_vector_free (b);

    gsl_matrix_free (Xe);
    gsl_matrix_free (A);
    gsl_matrix_free (XeX);
    gsl_matrix_free (Vaux);

    free (coefs);

    return(rc);
}

/**
 * @brief Clean up st_info
 *
 * @param st_info Pointer to container structure for Stata info
 * @return Frees memory allocated to st_info objects
 */
void ssf_free (struct StataInfo *st_info)
{
    if ( st_info->free >= 1 ) {
        free (st_info->groupvar);
        free (st_info->clustervar);
        free (st_info->X);
        free (st_info->y);
    }
    if ( st_info->free >= 2 ) {
        free (st_info->nobs);
        free (st_info->info);
    }
    if ( st_info->free >= 3 ) {
        free (st_info->cinfo);
    }
}
