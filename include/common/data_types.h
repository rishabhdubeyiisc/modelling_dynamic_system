#ifndef DATA_TYPES_H
#define DATA_TYPES_H

/**
 * @file data_types.h
 * @brief Core data structures for power system dynamics simulation
 * 
 * This header defines all the fundamental data structures used throughout
 * the power system simulation library including generator models, network
 * components, and state variables.
 */

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>

/* Note: These includes are needed here because this header defines structs using GSL types */

/* Legacy constant - use RVAL from common.h instead */
#define Rval 9

/**
 * @brief System constants and dimensions
 * 
 * Contains fundamental system dimensions used throughout the simulation.
 */
struct SystemConstants
{ 
    int NumberofGens;        /**< Number of generators in the system */
    int Number_of_lines;     /**< Number of transmission lines */
    int LOAD_FLOW_number;    /**< Number of buses in the system */
};
typedef struct SystemConstants SystemConstants;

/**
 * @brief Synchronous generator parameters
 * 
 * Contains all electrical and mechanical parameters for generator modeling.
 * Reactances: 0=steady-state, 1=transient, 2=sub-transient
 */
struct GeneratorParams
{   
    double  H;      /**< Inertia constant (seconds) */
    double  Xd0;    /**< d-axis synchronous reactance (pu) */
    double  Xd1;    /**< d-axis transient reactance (pu) */
    double  Xd2;    /**< d-axis sub-transient reactance (pu) */
    double  Xq0;    /**< q-axis synchronous reactance (pu) */
    double  Xq1;    /**< q-axis transient reactance (pu) */
    double  Xq2;    /**< q-axis sub-transient reactance (pu) */
    double  Tdo1;   /**< d-axis transient time constant (s) */
    double  Tdo2;   /**< d-axis sub-transient time constant (s) */
    double  Tqo1;   /**< q-axis transient time constant (s) */
    double  Tqo2;   /**< q-axis sub-transient time constant (s) */
};
typedef struct GeneratorParams GeneratorParams;

/**
 * @brief Exciter system parameters
 * 
 * Contains parameters for automatic voltage regulator (AVR) system.
 */
struct ExciterParams
{
    double  ka;     /**< Regulator gain */
    double  ta;     /**< Regulator time constant */
    double  ke;     /**< Exciter gain */
    double  te;     /**< Exciter time constant */
    double  kf;     /**< Stabilizing feedback gain */
    double  tf;     /**< Stabilizing feedback time constant */
};
typedef struct ExciterParams ExciterParams;


/**
 * @brief Transmission line parameters
 * 
 * Contains electrical parameters for transmission line modeling.
 */ 
struct TransmissionLineParams
{   
    int     bus1;   /**< From bus number */
    int     bus2;   /**< To bus number */
    double  R;      /**< Resistance (pu) */
    double  X;      /**< Reactance (pu) */
    double  B;      /**< Susceptance (pu) */
};
typedef struct TransmissionLineParams TransmissionLineParams;


/**
 * @brief Bus data from load flow analysis
 * 
 * Contains power flow results and bus specifications.
 */
struct BusData
{   
    int     bus_number; /**< Bus number identifier */
    int     bus_type;   /**< Bus type (1=slack, 2=PV, 3=PQ) */
    double  V;          /**< Voltage magnitude (pu) */
    double  theta;      /**< Voltage angle (degrees) */
    double  Pg;         /**< Generated active power (pu) */
    double  Qg;         /**< Generated reactive power (pu) */
    double  Pl;         /**< Load active power (pu) */
    double  Ql;         /**< Load reactive power (pu) */
};
typedef struct BusData BusData;

/**
 * @brief Complete network data structure
 * 
 * Contains all system data including generators, lines, and buses.
 */
struct NetworkData
{   
    SystemConstants         * constants;        /**< System dimensions */
    GeneratorParams         * Generator_ps;     /**< Generator parameters */
    ExciterParams           * exciter_ps;       /**< Exciter parameters */
    BusData                 * Load_flow_ps;     /**< Bus data */
    TransmissionLineParams  * trans_line_para;  /**< Transmission line data */
};
typedef struct NetworkData NetworkData;


/**
 * @brief Initial conditions for generator states
 * 
 * Contains steady-state initial values for all generator variables.
 */
struct InitialConditions
{ 
    gsl_complex I_0;            /**< Initial generator current */
    gsl_complex Eq_0;           /**< Initial EMF behind transient reactance */
    gsl_complex id_0;           /**< Initial d-axis current */
    gsl_complex iq_0;           /**< Initial q-axis current */
    gsl_complex vd_0;           /**< Initial d-axis voltage */
    gsl_complex vq_0;           /**< Initial q-axis voltage */
    gsl_complex Efd_0;          /**< Initial field voltage */
    gsl_complex Eq_dash_0;      /**< Initial transient EMF q-axis */
    gsl_complex Ed_dash_0;      /**< Initial transient EMF d-axis */
    gsl_complex Pe_0;           /**< Initial electrical power */
    gsl_complex Pm_0;           /**< Initial mechanical power */
    
    double  Ed_das_das;         /**< Initial sub-transient EMF d-axis */
    double  Eq_das_das;         /**< Initial sub-transient EMF q-axis */
    double  delta_0;            /**< Initial rotor angle */
    double  VQ;                 /**< Initial Q-axis bus voltage */
    double  VD;                 /**< Initial D-axis bus voltage */
    
    gsl_complex I_0_DQFRAME;    /**< Initial current in DQ frame */
};
typedef struct InitialConditions InitialConditions;

/**
 * @brief Network admittance matrix
 * 
 * Contains complex admittance matrix for network analysis.
 */
struct AdmittanceMatrix
{
    gsl_complex MAT[Rval][Rval];  /**< Complex admittance matrix */
};
typedef struct AdmittanceMatrix AdmittanceMatrix;


/**
 * @brief Generator dynamic state variables
 * 
 * Contains all time-varying state variables for generator dynamics.
 */
struct GeneratorStates
{
    double  Ed_das;         /**< Transient EMF d-axis */
    double  Eq_das;         /**< Transient EMF q-axis */
    double  Ed_das_das;     /**< Sub-transient EMF d-axis */
    double  Eq_das_das;     /**< Sub-transient EMF q-axis */
    double  delta;          /**< Rotor angle (rad) */
    double  slip;           /**< Rotor speed deviation (pu) */
    double  Efd;            /**< Field voltage */
    double  E_dummy;        /**< Auxiliary EMF state */
};
typedef struct GeneratorStates GeneratorStates;

/**
 * @brief Network bus voltage states
 * 
 * Contains bus voltages in DQ reference frame.
 */
struct NetworkStates
{   
    double  VQ;     /**< Q-axis bus voltage */
    double  VD;     /**< D-axis bus voltage */
};
typedef struct NetworkStates NetworkStates;

/**
 * @brief Current in DQ coordinates
 * 
 * Contains generator current components in DQ reference frame.
 */
struct CurrentDQ
{
    double  IQ;     /**< Q-axis current component */
    double  ID;     /**< D-axis current component */
};
typedef struct CurrentDQ CurrentDQ;


/**
 * @brief Time derivatives of generator state variables
 * 
 * Contains the time derivatives for numerical integration.
 */
struct StateDerivatives 
{
    double  f_of_Ed_das;        /**< dEd'/dt */
    double  f_of_Eq_das;        /**< dEq'/dt */
    double  f_of_Ed_das_das;    /**< dEd''/dt */
    double  f_of_Eq_das_das;    /**< dEq''/dt */
    double  f_of_delta;         /**< dδ/dt */
    double  f_of_slip;          /**< ds/dt */
    double  f_of_Efd;           /**< dEfd/dt */
    double  f_of_E_dummy;       /**< dE_dummy/dt */
};
typedef struct StateDerivatives StateDerivatives;

struct Y_block_2_2
{
    double Y_2X2[2][2];
};
typedef struct Y_block_2_2 Y_2X2;


#endif



