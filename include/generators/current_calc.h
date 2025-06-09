#ifndef CURRENT_CALC_H
#define CURRENT_CALC_H

#include "../common/types_fwd.h"

/**
 * Update generator currents based on states and terminal voltages
 * @param X_iter Generator state variables
 * @param vd_iter d-axis terminal voltage
 * @param vq_iter q-axis terminal voltage
 * @param ptr_id Pointer to d-axis current output
 * @param ptr_iq Pointer to q-axis current output
 * @param Xd2_iter d-axis sub-transient reactance
 */
void current_update_new(GeneratorStates X_iter,
                        double vd_iter,
                        double vq_iter,
                        double * ptr_id,
                        double * ptr_iq,
                        double Xd2_iter);

#endif 