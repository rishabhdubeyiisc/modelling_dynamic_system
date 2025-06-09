#ifndef POWER_SYSTEM_API_H
#define POWER_SYSTEM_API_H

/**
 * @file power_system.h
 * @brief Master header for Power System Dynamics Simulation Library
 * @version 1.0
 * @date 2025
 * 
 * This is the main header file for the power system dynamics simulation library.
 * Include this single header to access all public APIs for:
 * - Multi-machine synchronous generator modeling
 * - Network admittance matrix operations
 * - Fault simulation and analysis
 * - Time-domain dynamics simulation
 * - Mathematical utilities and linear algebra
 * 
 * @author Power Systems Research Group
 * 
 * @example
 * ```c
 * #include "api/power_system.h"
 * 
 * int main() {
 *     // Read system data
 *     NetworkData data = loadSystemData();
 *     
 *     // Initialize generators
 *     InitialConditions* init_states = initializeGenerator(data);
 *     
 *     // Run simulation
 *     partitioned_solver_sm(data, init_states, 0.002, 600, ...);
 *     
 *     return 0;
 * }
 * ```
 */

/* Core system includes */
#include "../common/common.h"

/* Public API modules */
#include "simulation.h"    /* Main simulation interface */
#include "network.h"       /* Network analysis functions */
#include "generators.h"    /* Generator modeling functions */
#include "math.h"          /* Mathematical utilities */

/* User interface modules */
#include "../io/user_input.h"      /* Parameter input handling */
#include "../core/fault_config.h"  /* Fault configuration interface */

/**
 * @brief Library version information
 */
#define POWER_SYSTEM_LIB_VERSION_MAJOR 1
#define POWER_SYSTEM_LIB_VERSION_MINOR 0
#define POWER_SYSTEM_LIB_VERSION_PATCH 0

/**
 * @brief Get library version string
 * @return Version string in format "major.minor.patch"
 */
static inline const char* power_system_version(void) {
    return "1.0.0";
}

/**
 * @brief Library feature flags
 */
#define POWER_SYSTEM_HAS_FAULT_SIMULATION 1
#define POWER_SYSTEM_HAS_MULTI_MACHINE 1
#define POWER_SYSTEM_HAS_INTERACTIVE_UI 1
#define POWER_SYSTEM_HAS_CSV_OUTPUT 1

#endif /* POWER_SYSTEM_API_H */ 