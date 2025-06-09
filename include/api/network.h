#ifndef NETWORK_API_H
#define NETWORK_API_H

#include "../common/common.h"

/**
 * @file network.h
 * @brief Public API for power system network analysis
 * 
 * This header provides interfaces for network matrix operations,
 * fault analysis, and power flow calculations.
 */

/**
 * @brief Convert complex Y-bus matrix to real 2x2 block format
 * @param Y Complex Y-bus matrix
 * @param All_data Network data structure for sizing
 * @return Dynamically allocated real matrix in 2x2 block format
 * 
 * Converts complex admittance matrix Y = G + jB into real matrix format:
 * [G -B]
 * [B  G]
 * This transformation enables real-valued network equation solving.
 */
double** convertComplexToReal(AdmittanceMatrix Y, NetworkData All_data);

/**
 * @brief Create Y-bus matrix with line fault applied
 * @param line_set Line index for fault (1-based indexing)
 * @param Data Network data structure
 * @return Y-bus matrix with specified line fault
 * 
 * Modifies the network admittance matrix to represent a fault condition
 * on the specified transmission line by splitting the line impedance.
 */
AdmittanceMatrix createFaultMatrix(int line_set, NetworkData Data);

/**
 * @brief Solve network equations for given injections
 * @param DATA Network data structure
 * @param Number_of_generators Number of generators in system
 * @param Number_of_buses Number of buses in system
 * @param pointer_nw_vector Output: network voltage states (VQ, VD)
 * @param ptr_X_vector Generator state variables
 * @param Z_AUG Network impedance matrix (inverse of Y_AUG)
 * @param ptr_iq q-axis generator currents
 * @param ptr_id d-axis generator currents  
 * @param ptr_Vt Terminal voltage magnitudes (output)
 * 
 * Solves the network equations V = Z * I where I includes generator
 * Norton equivalent currents and load currents. Updates network
 * voltage states and terminal voltages.
 */
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

/**
 * @brief Print complex matrix for debugging
 * @param Y Complex matrix to print
 * @param size Matrix dimension
 * 
 * Utility function for debugging complex matrix values.
 */
void print_matrix_complex(AdmittanceMatrix Y, int size);

#endif /* NETWORK_API_H */ 