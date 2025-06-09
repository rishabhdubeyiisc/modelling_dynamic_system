#ifndef NETWORK_FWD_H
#define NETWORK_FWD_H

/**
 * @file network_fwd.h
 * @brief Fast-compiling network interface using forward declarations
 * 
 * This header provides network function declarations without heavy includes.
 * Use this when you only need to declare functions, not implement them.
 */

#include "../common/types_fwd.h"

/* Network analysis functions - forward declarations only */
double** convertComplexToReal(AdmittanceMatrix Y, NetworkData All_data);
AdmittanceMatrix createFaultMatrix(int line_set, NetworkData Data);

void solveNetworkStep(
    NetworkData DATA,
    int Number_of_generators, 
    int Number_of_buses,
    NetworkStates* pointer_nw_vector, 
    GeneratorStates* ptr_X_vector,
    double** Z_AUG,
    double* ptr_iq,
    double* ptr_id,
    double* ptr_Vt
);

void print_matrix_complex(AdmittanceMatrix Y, int size);

#endif /* NETWORK_FWD_H */ 