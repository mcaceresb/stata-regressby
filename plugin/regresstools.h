#ifndef REGRESSTOOLS_H
#define REGRESSTOOLS_H    1

#define BUF_MAX 4096

void sf_printf (const char *fmt, ...);
void sf_running_timer (clock_t *timer, const char *msg);

void sf_gsl_printf_vector(const char * fmt, const gsl_vector * v);
void sf_gsl_printf_matrix(const char * fmt, const gsl_matrix * M);

ST_retcode sf_scalar_size (char *st_scalar, GT_size *sval);
ST_retcode sf_oom_error (char *step_desc, char *obj_desc);

ST_retcode sf_gsl_get_variable(
    gsl_vector * v,
    GT_size k,
    GT_size in1,
    GT_size in2
);

ST_retcode sf_gsl_get_varlist(
    gsl_matrix * M,
    GT_size k1,
    GT_size k2,
    GT_size in1,
    GT_size in2
);

GT_size sf_panelsetup (GT_size *group, GT_size *info, GT_size N);

#endif /* regresstools.h */
