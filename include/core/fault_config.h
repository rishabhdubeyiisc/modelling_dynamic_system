#ifndef FAULT_CONFIG_H
#define FAULT_CONFIG_H

#include "../common/types_fwd.h"

/**
 * Configure fault simulation parameters interactively
 * @param Net_data_main Network data structure
 * @param Y_complex_aug_main Healthy Y-bus matrix
 * @param Z_AUG_fault Pointer to fault impedance matrix (will be allocated)
 * @param fault_start Pointer to fault start time
 * @param fault_end Pointer to fault end time
 * @return 1 if fault enabled, 0 otherwise
 */
int configure_fault_simulation(NetworkData Net_data_main, AdmittanceMatrix Y_complex_aug_main,
                              double ***Z_AUG_fault, double *fault_start, double *fault_end);

#endif 