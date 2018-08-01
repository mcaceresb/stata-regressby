#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "gttypes.h"
#include "stplugin.h"
#include "regresstools.h"

/**
 * @brief Short wrapper to print to Stata
 *
 * Basic wrapper to print formatted strings to Stata
 *
 * @param *fmt a string to format
 * @param ... Arguments to pass to pritnf
 * @return Prints to Stata's console
 */
void sf_printf (const char *fmt, ...)
{
    va_list args;
    va_start (args, fmt);
    char buf[BUF_MAX];
    vsprintf (buf, fmt, args);
    SF_display (buf);
    // printf (buf);
    va_end (args);
}

/**
 * @brief Short wrapper to print error to Stata
 *
 * Basic wrapper to print formatted error strings to Stata
 *
 * @param *fmt a string to format
 * @param ... Arguments to pass to pritnf
 * @return Prints to Stata's console
 */
void sf_errprintf (const char *fmt, ...)
{
    va_list args;
    va_start (args, fmt);
    char buf[BUF_MAX];
    vsprintf (buf, fmt, args);
    SF_error (buf);
    va_end (args);
}

/**
 * @brief Update a running timer and print a message to satata console
 *
 * Prints a messasge to Stata that the running timer @timer was last set
 * @diff seconds ago. It then updates the timer to the current time.
 *
 * @param timer clock object containing time since last udpate
 * @param msg message to print before # of seconds
 * @return Print time since last update to Stata console
 */
void sf_running_timer (clock_t *timer, const char *msg)
{
    double diff  = (double) (clock() - *timer) / CLOCKS_PER_SEC;
    sf_printf (msg);
    sf_printf ("; %.3f seconds.\n", diff);
    *timer = clock();
}

/**
 * @brief Read scalar into unsigned integer
 *
 * @param st_scalar name of Stata scalar
 * @param sval Scalar value
 * @return Read scalar into GT_size variable
 */
ST_retcode sf_scalar_size (char *st_scalar, GT_size *sval)
{
    ST_double _double;
    ST_retcode rc = SF_scal_use(st_scalar, &_double);
    if ( rc == 0 ) {
        *sval = (GT_size) _double;
    }
    return (rc);
}

/**
 * @brief Panel setup based on vector
 */
ST_retcode sf_oom_error (char *step_desc, char *obj_desc)
{
    sf_errprintf ("%s: Unable to allocate memory for object '%s'.\n", step_desc, obj_desc);
    return (1702);
}

/**
 * @brief Parse variable into vector from Stata
 *
 * Loop through a variable's observations and parse it into a GSL
 * vector. @v must be the length of the number of observations to be
 * read.
 *
 * @param v GSL vector to save into.
 * @param k Position of variable in varlist passed to C
 * @param in1 First observation to be read
 * @param in2 Last observation to be read
 * @return Copies @k variable into @v
 *
 * @see sf_gsl_get_varlist
 */
ST_retcode sf_gsl_get_variable(
    gsl_vector * v,
    GT_size k,
    GT_size in1,
    GT_size in2)
{
    ST_retcode rc ;
    ST_double  z ;
    GT_size i;
    for (i = in1; i <= in2; i++) {
        if ( (rc = SF_vdata(k, i, &z)) ) return(rc);
        gsl_vector_set (v, i - in1, z);
    }
    return(0);
}

/**
 * @brief Parse varlist into matrix from Stata
 *
 * Loop through a varlist's observations and parse it into a GSL
 * matrix. @M must be the length of the number of observations to be
 * read and width of the number of variables to be read.
 *
 * @param M GSL matrix to save into.
 * @param k1 Position of variable in varlist passed to C
 * @param k2 Position of variable in varlist passed to C
 * @param in1 First observation to be read
 * @param in2 Last observation to be read
 * @return Copies variables @k1 to @k2 into @M
 *
 * @see sf_gsl_get_varlist
 */
ST_retcode sf_gsl_get_varlist(
    gsl_matrix * M,
    GT_size k1,
    GT_size k2,
    GT_size in1,
    GT_size in2)
{
    ST_retcode rc ;
    ST_double  z ;
    GT_size i, j;
    for (i = in1; i <= in2; i++) {
        for (j = k1; j <= k2; j++) {
            if ( (rc = SF_vdata(j, i, &z)) ) return(rc);
            gsl_matrix_set (M, i - in1, j - k1, z);
        }
    }
    return(0);
}

/**
 * @brief Set up variables for panel hashes
 *
 * Using a sorted integer variable, generate info array with start and
 * ending positions of each group in the sorted vector.
 *
 * @param group Array of sorted integers
 * @info Where to store starting and ending positions
 * @J Number of groups
 *
 * @return info arary with start and end positions of each group
 */
GT_size sf_panelsetup (GT_size *group, GT_size *info, GT_size N)
{
    GT_size *gcmp, *gptr, i, j;

    j = 0;
    i = 0;
    info[j] = 0;
    gcmp    = group;

    for (gptr = group + 1; gptr < group + N; gptr++, i++) {
        if ( *gcmp != *gptr ) {
            j++;
            *gcmp   = *gptr;
            info[j] = i + 1;
        }
    }

    j++;
    info[j] = N;

    return (j);
} 

/**
 * @brief Print gsl vector into Stata
 *
 * @param v gsl vector to print
 * @param fmt printf format
 *
 * @see sf_gsl_printf_matrix
 */
void sf_gsl_printf_vector(const char * fmt, const gsl_vector * v)
{
    GT_size i;
    for (i = 0; i < v->size; i++) {
        sf_printf(fmt, gsl_vector_get(v, i));
    }
}

/**
 * @brief Print gsl matrix into Stata (rows are ended with "\n")
 *
 * @param M gsl matrix to print
 * @param fmt printf format
 *
 * @see sf_gsl_printf_vector
 */
void sf_gsl_printf_matrix(const char * fmt, const gsl_matrix * M)
{
    GT_size i, j;
    for (i = 0; i < M->size1; i++) {
        for (j = 0; j < M->size2; j++) {
            sf_printf(fmt, gsl_matrix_get(M, i, j));
        }
        sf_printf("\n");
    }
}
