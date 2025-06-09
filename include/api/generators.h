#ifndef GENERATORS_API_H
#define GENERATORS_API_H

#include "../common/common.h"

/**
 * @file generators.h
 * @brief Public API for synchronous generator modeling
 * 
 * This header provides interfaces for generator electrical calculations,
 * power computations, and current updates in power system dynamics.
 */

/**
 * @brief Calculate electrical power output of synchronous generator
 * @param Ed2 d-axis sub-transient EMF (pu)
 * @param id d-axis stator current (pu)
 * @param Eq2 q-axis sub-transient EMF (pu)  
 * @param iq q-axis stator current (pu)
 * @param Xd2 d-axis sub-transient reactance (pu)
 * @param Xq2 q-axis sub-transient reactance (pu)
 * @return Electrical power output in per unit
 * 
 * Computes electrical power using the Park's transformation equations:
 * Pe = Ed2*id + Eq2*iq + (Xd2-Xq2)*id*iq
 */
double computeElectricalPower(double Ed2, double id, double Eq2, double iq, double Xd2, double Xq2);

/**
 * @brief Update generator currents from terminal voltages and states
 * @param X_iter Current generator state variables
 * @param vd_iter d-axis terminal voltage (pu)
 * @param vq_iter q-axis terminal voltage (pu)
 * @param ptr_id Output: d-axis current (pu)
 * @param ptr_iq Output: q-axis current (pu)
 * @param Xd2_iter d-axis sub-transient reactance (pu)
 * 
 * Updates generator stator currents based on terminal voltages and
 * internal EMFs using generator equivalent circuit equations:
 * I = (E - V) / jXd2
 */
void current_update_new(
    GeneratorStates X_iter,
    double vd_iter,
    double vq_iter,
    double* ptr_id,
    double* ptr_iq,
    double Xd2_iter
);

/**
 * @brief Calculate mechanical torque from power and slip
 * @param power Electrical power (pu)
 * @param slip Generator slip (pu)
 * @return Mechanical torque in per unit
 * 
 * Helper function for generator mechanical dynamics calculations.
 */
double Torque_calc(double power, double slip);

/**
 * @brief Euler forward integration for generator states (without exciter)
 * @param current_state Current generator states
 * @param next_state Output: next time step states
 * @param derivatives State derivatives
 * @param del_t Time step size (seconds)
 * 
 * Performs numerical integration of generator differential equations
 * using Euler forward method for cases without exciter dynamics.
 */
void Euler_forward_wo_Exiter(
    GeneratorStates* current_state,
    GeneratorStates* next_state, 
    StateDerivatives* derivatives,
    double del_t
);

/**
 * @brief Calculate state derivatives for single generator
 * @param gen_index Generator index
 * @param derivatives Output: computed derivatives
 * @param All_data Network data structure
 * @param states Current generator states
 * @param vq q-axis terminal voltage
 * @param vd d-axis terminal voltage
 * @param iq q-axis current
 * @param id d-axis current
 * @param del_t Time step size
 * 
 * Computes time derivatives of generator state variables including
 * rotor angle, speed, and EMF dynamics for numerical integration.
 */
void F_VECTOR_CALC_single(
    int gen_index,
    StateDerivatives* derivatives,
    NetworkData All_data,
    GeneratorStates states,
    double vq,
    double vd,
    double iq,
    double id,
    double del_t
);

#endif /* GENERATORS_API_H */ 