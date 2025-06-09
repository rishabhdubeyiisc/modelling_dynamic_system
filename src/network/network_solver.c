#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "../../include/network/network_solver.h"
#include "../../include/common/data_types.h"

//WORKS_ON_NETWROK_VECTOR_ONLY
void solveNetworkStep(NetworkData DATA,
                      int Number_of_generators, 
                      int Number_of_buses,
                      NetworkStates * pointer_nw_vector, 
                      GeneratorStates * ptr_X_vector,
                      double ** Z_AUG,
                      double * ptr_iq,
                      double * ptr_id,
                      double * ptr_Vt) {       

    double Xd2_local[Number_of_generators];
    double Xq2_local[Number_of_generators];

    double i_Q_network[Number_of_generators];
    double i_D_network[Number_of_generators];
    
    double I_inject_NW      [2*Number_of_buses];
    double V_vec            [2*Number_of_buses];

    double iq_loc [Number_of_generators];
   
    // ✅ FIX: Use local indices instead of modifying input pointers
    for (int i = 0; i < Number_of_generators; i++) {
        Xd2_local[i]=DATA.Generator_ps[i].Xd2;
        Xq2_local[i]=DATA.Generator_ps[i].Xq2;
        iq_loc[i]=ptr_iq[i]; // ✅ Use array indexing instead of pointer increment
    }

    // ✅ FIX: Use local index instead of modifying input pointer
    for (int i = 0; i < Number_of_generators; i++) {   
        double delta_L=ptr_X_vector[i].delta; // ✅ Use array indexing

        ptr_X_vector[i].E_dummy=(-1)*(Xq2_local[i]-Xd2_local[i])*(iq_loc[i]);

        double Ed2_add_Edummy=ptr_X_vector[i].E_dummy + ptr_X_vector[i].Ed_das_das;
        
        gsl_complex numerator = gsl_complex_rect(ptr_X_vector[i].Eq_das_das,Ed2_add_Edummy);
        
        gsl_complex Gen_impedance = gsl_complex_rect(0,Xd2_local[i]);
        
        gsl_complex I_gen_norton_dq = gsl_complex_div(numerator,Gen_impedance);

        // //calculatr notn in ab
        gsl_complex delta_complex = gsl_complex_rect(0,delta_L);

        gsl_complex shifter = gsl_complex_exp(delta_complex);

        gsl_complex I_norton_shifted = gsl_complex_mul(I_gen_norton_dq,shifter);

        i_Q_network[i]=(I_norton_shifted.dat[0]);

        i_D_network[i]=(I_norton_shifted.dat[1]);
        //initialize the V gens to zero 1st 

        I_inject_NW[2*i]        =i_Q_network[i];
        I_inject_NW[2*i+1]      =i_D_network[i];
    }

    // ✅ FIX: Use proper generalized bounds for load buses
    for (int i = 2*Number_of_generators; i < 2*Number_of_buses; i++) {
        I_inject_NW[i]=0;
    }
    
    //I_vec_lenght 18
    double sum = 0;
    for (int i = 0; i < 2*Number_of_buses; i++) {       
        for (int j = 0; j < 2*Number_of_buses; j++) {
            double prod = 0;
            prod = Z_AUG[i][j]*I_inject_NW[j];
            sum = sum + prod ; 
        }
        V_vec[i]=sum;
        sum = 0;
    }
    
    // ✅ FIX: Use local index instead of modifying input pointer
    for (int i = 0; i < Number_of_buses; i++) {
          pointer_nw_vector[i].VQ = V_vec[2*i];
          pointer_nw_vector[i].VD = V_vec[2*i + 1];  
    }
} 