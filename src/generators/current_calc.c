#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "../../include/generators/current_calc.h"
#include "../../include/common/data_types.h"

void current_update_new(GeneratorStates X_iter,
                        double vd_iter,
                        double vq_iter,
                        double * ptr_id,
                        double * ptr_iq,
                        double Xd2_iter) {
    double Eq2 = X_iter.Eq_das_das;
    double Ed2 = X_iter.Ed_das_das;
    double E_dummy = X_iter.E_dummy;
    double Ed_add_Edum = Ed2 + E_dummy;
    double vq = vq_iter;
    double vd = vd_iter;
    
    /* Numerical robustness: ensure minimum reactance to avoid division by zero */
    if (fabs(Xd2_iter) < 1e-8) {
        Xd2_iter = (Xd2_iter >= 0) ? 1e-8 : -1e-8;
    }
    
    gsl_complex v_park = gsl_complex_rect(vq, vd);
    gsl_complex Z = gsl_complex_rect(0, Xd2_iter);

    gsl_complex Eq_Ed_Edum = gsl_complex_rect(Eq2, Ed_add_Edum);

    gsl_complex numerator = gsl_complex_sub(Eq_Ed_Edum, v_park);

    gsl_complex I = gsl_complex_div(numerator, Z);

    *ptr_iq = I.dat[0];
    *ptr_id = I.dat[1];
    
    /* Numerical robustness: limit extreme current values */
    if (fabs(*ptr_iq) > 10.0) {
        *ptr_iq = (*ptr_iq > 0) ? 10.0 : -10.0;
    }
    if (fabs(*ptr_id) > 10.0) {
        *ptr_id = (*ptr_id > 0) ? 10.0 : -10.0;
    }
} 