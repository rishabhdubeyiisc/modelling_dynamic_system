#ifndef COMMON_H
#define COMMON_H

/* Standard C library includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* GSL library includes */
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

/* Core data types */
#include "data_types.h"

/* System-wide constants */
#define PI acos(-1.0)
#define OMEGA_BASE (2*PI*50)        /* Base angular frequency (rad/s) */
#define T_dummy 0.1                 /* Dummy time constant */
#define ONE_BY_T_DUMMY 10           /* Inverse of dummy time constant */
#define ERROR 0.0001                /* Numerical error tolerance */
#define ELR_FACTOR 0.5              /* Error learning rate factor */

/* System limits */
#ifndef NG_MAX
#define NG_MAX 20                   /* Maximum number of generators */
#endif

#ifndef RVAL  
#define RVAL 9                      /* Matrix dimension constant */
#endif

/* Utility macros */
#define DEG_TO_RAD(deg) ((deg) * PI / 180.0)
#define RAD_TO_DEG(rad) ((rad) * 180.0 / PI)

#endif /* COMMON_H */ 