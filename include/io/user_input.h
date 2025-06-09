#ifndef USER_INPUT_H
#define USER_INPUT_H

#include "../common/types_fwd.h"

/**
 * Get simulation parameters from user input
 * @param target_g Pointer to target generator index
 * @param vref_step_delta Pointer to Vref step change value
 * @param pm_step_delta Pointer to Pm step change value  
 * @param vref_step_time Pointer to Vref step time
 * @param pm_step_time Pointer to Pm step time
 * @param Initial_state_main Initial state data for generators
 */
void get_simulation_parameters(int *target_g, double *vref_step_delta, 
                              double *pm_step_delta, double *vref_step_time, 
                              double *pm_step_time, InitialConditions *Initial_state_main);

#endif 