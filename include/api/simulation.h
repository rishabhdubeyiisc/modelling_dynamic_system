#ifndef SIMULATION_API_H
#define SIMULATION_API_H

#include "../common/common.h"

/**
 * @file simulation.h
 * @brief Public API for power system dynamics simulation
 * 
 * This header provides the main simulation interfaces for power system
 * dynamics analysis including data reading, initialization, and solver functions.
 */

/**
 * @brief Read system data from files
 * @return Combined structure containing all network data
 * 
 * Reads generator parameters, load flow data, and transmission line
 * parameters from the data directory.
 */
NetworkData loadSystemData(void);

/**
 * @brief Calculate initial generator states from load flow data
 * @param NetData Network data structure containing load flow and generator data
 * @return Array of initial generator states
 * 
 * Computes steady-state initial conditions for all generators including
 * EMFs, currents, and state variables.
 */
InitialConditions* initializeGenerator(NetworkData NetData);

/**
 * @brief Create augmented Y-bus matrix for the system
 * @param Data Network data structure
 * @return Y-bus matrix in complex format
 * 
 * Constructs the network admittance matrix including generator and
 * transmission line contributions.
 */
AdmittanceMatrix createAugmentedAdmittanceMatrix(NetworkData Data);

/**
 * @brief Main partitioned solver for multi-machine dynamics
 * @param All_data Network data structure
 * @param Initial_state Array of initial generator states
 * @param del_t Time step size (seconds)
 * @param END_TIME Simulation end time (seconds)
 * @param Z_AUG_healthy Healthy system impedance matrix
 * @param Z_AUG_fault Fault system impedance matrix (NULL if no fault)
 * @param fault_enabled 1 if fault simulation enabled, 0 otherwise
 * @param fault_start Fault start time (seconds)
 * @param fault_end Fault end time (seconds)
 * @param TARGET_G Target generator index for control actions
 * @param VREF_STEP_TIME Time of Vref step change (seconds)
 * @param VREF_STEP_DELTA Magnitude of Vref step change (pu)
 * @param PM_STEP_TIME Time of mechanical power step change (seconds)
 * @param PM_STEP_DELTA Magnitude of mechanical power step change (pu)
 * 
 * Performs time-domain simulation of multi-machine power system dynamics
 * with optional fault simulation and control input steps.
 */
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

#endif /* SIMULATION_API_H */ 