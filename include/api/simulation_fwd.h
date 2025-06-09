#ifndef SIMULATION_FWD_H
#define SIMULATION_FWD_H

/**
 * @file simulation_fwd.h
 * @brief Fast-compiling simulation interface using forward declarations
 * 
 * This header provides simulation function declarations without heavy includes.
 * Use this when you only need to declare functions, not implement them.
 * For implementation, include the full "api/simulation.h" instead.
 */

#include "../common/types_fwd.h"

/* Core simulation functions - forward declarations only */
NetworkData loadSystemData(void);
InitialConditions* initializeGenerator(NetworkData NetData);
AdmittanceMatrix createAugmentedAdmittanceMatrix(NetworkData Data);

void partitioned_solver_sm(
    NetworkData All_data,
    InitialConditions* Initial_state,
    double del_t,
    double END_TIME,
    double** Z_AUG_healthy,
    double** Z_AUG_fault,
    int fault_enabled,
    double fault_start,
    double fault_end,
    int TARGET_G,
    double VREF_STEP_TIME,
    double VREF_STEP_DELTA,
    double PM_STEP_TIME,
    double PM_STEP_DELTA
);

#endif /* SIMULATION_FWD_H */ 