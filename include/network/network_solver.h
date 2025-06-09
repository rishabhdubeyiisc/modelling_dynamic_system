#ifndef NETWORK_SOLVER_H
#define NETWORK_SOLVER_H

#include "../common/types_fwd.h"

/**
 * Network solver for power system dynamics simulation
 * @param DATA Network data structure
 * @param Number_of_generators Number of generators in the system
 * @param Number_of_buses Number of buses in the system  
 * @param pointer_nw_vector Network state variables output
 * @param ptr_X_vector Generator state variables
 * @param Z_AUG Network impedance matrix
 * @param ptr_iq q-axis current array
 * @param ptr_id d-axis current array
 * @param ptr_Vt Terminal voltage array
 */
void solveNetworkStep(NetworkData DATA,
                      int Number_of_generators, 
                      int Number_of_buses,
                      NetworkStates * pointer_nw_vector, 
                      GeneratorStates * ptr_X_vector,
                      double ** Z_AUG,
                      double * ptr_iq,
                      double * ptr_id,
                      double * ptr_Vt);

#endif 