#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "../common/data_types.h"

/**
 * @brief Partitioned time-domain solver (production)
 * 
 * Solves generator dynamics and network equations in a partitioned manner:
 * 1. Solve network algebraically at each time step
 * 2. Integrate generator dynamics using RK4
 * 
 * Fast and robust for most power system studies.
 */
void partitioned_solver_sm(
        NetworkData All_data,
        InitialConditions   *Initial_state,
        double     del_t,
        double     END_TIME,
        double   **Z_AUG_healthy,
        double   **Z_AUG_fault,
        int        fault_enabled,
        double     fault_start,
        double     fault_end,
        int        TARGET_G,
        double     VREF_STEP_TIME,
        double     VREF_STEP_DELTA,
        double     PM_STEP_TIME,
        double     PM_STEP_DELTA);

/**
 * @brief Jacobian-based simultaneous solver (experimental)
 * 
 * Solves generator and network equations simultaneously using:
 * 1. Implicit Euler integration
 * 2. Newton-Raphson with full system Jacobian
 * 
 * More accurate for tightly-coupled systems but slower.
 */
void jacobian_solver_sm(
        NetworkData All_data,
        InitialConditions   *Initial_state,
        double     del_t,
        double     END_TIME,
        double   **Z_AUG_healthy,
        double   **Z_AUG_fault,
        int        fault_enabled,
        double     fault_start,
        double     fault_end,
        int        TARGET_G,
        double     VREF_STEP_TIME,
        double     VREF_STEP_DELTA,
        double     PM_STEP_TIME,
        double     PM_STEP_DELTA);

#endif /* SIMULATOR_H */ 