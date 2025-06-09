#ifndef TYPES_FWD_H
#define TYPES_FWD_H

/**
 * @file types_fwd.h
 * @brief Forward declarations for core data types
 * 
 * This header provides forward declarations without full type definitions,
 * allowing headers to declare function signatures without including the
 * heavy data_types.h file. This reduces compilation dependencies.
 */

/* Forward declarations - use these when you only need to declare pointers or references */
typedef struct SystemConstants SystemConstants;
typedef struct GeneratorParams GeneratorParams;
typedef struct ExciterParams ExciterParams;
typedef struct TransmissionLineParams TransmissionLineParams;
typedef struct BusData BusData;
typedef struct NetworkData NetworkData;
typedef struct InitialConditions InitialConditions;
typedef struct AdmittanceMatrix AdmittanceMatrix;
typedef struct GeneratorStates GeneratorStates;
typedef struct NetworkStates NetworkStates;
typedef struct CurrentDQ CurrentDQ;
typedef struct StateDerivatives StateDerivatives;
typedef struct Y_block_2_2 Y_2X2;

/* Note: GSL types (gsl_complex, etc.) must be included from GSL headers when needed */

#endif /* TYPES_FWD_H */ 