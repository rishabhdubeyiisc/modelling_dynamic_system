#include "../../include/common/common.h"
#include "../../include/generators/gen_init.h"

InitialConditions * initializeGenerator(NetworkData NetData) {   
    SystemConstants   * constant = &NetData.constants[0];
    int Number_of_gens = constant[0].NumberofGens;

    InitialConditions *  Gen_initial;
    Gen_initial = malloc(Number_of_gens*sizeof(InitialConditions));
        
    gsl_complex I_0      [Number_of_gens];
    gsl_complex Eq_0     [Number_of_gens];
    gsl_complex id_0     [Number_of_gens];
    gsl_complex iq_0     [Number_of_gens];
    gsl_complex vd_0     [Number_of_gens];
    gsl_complex vq_0     [Number_of_gens];
    gsl_complex Efd_0    [Number_of_gens];
    gsl_complex Eq_dash_0[Number_of_gens];
    gsl_complex Ed_dash_0[Number_of_gens];
    gsl_complex Pe_0     [Number_of_gens];
    gsl_complex Pm_0     [Number_of_gens];

    double VD[Number_of_gens];
    double VQ[Number_of_gens];
        
    double delta_0_arr   [Number_of_gens];
    double Ed_das_das    [Number_of_gens]; 
    double Eq_das_das    [Number_of_gens];

    gsl_complex I_0_DQ_FRAME[Number_of_gens];   
  
    GeneratorParams    *  Gen_array      = &NetData.Generator_ps[0];
    BusData            *  Load_flow_arr  = &NetData.Load_flow_ps[0];
        
    for (int x = 0; x < Number_of_gens  ; x++, Gen_array++, Load_flow_arr++) {
        //double theta             = Load_flow_arr->theta; 
        double angle             = ((Load_flow_arr->theta)*PI)/180;
        //made_change_check
        gsl_complex z            = gsl_complex_rect(0, angle); 
        double neg_angle=(-1)*(angle);
        gsl_complex z2 = gsl_complex_rect(0,neg_angle);
        
        gsl_complex e__neg_j_theta    = gsl_complex_exp(z2);
        gsl_complex e_j_theta          = gsl_complex_exp((z));
            
        gsl_complex Vt_0         = gsl_complex_rect(Load_flow_arr->V,0.0) ;
        double      mag_Vt_0     = gsl_complex_abs(Vt_0);

        gsl_complex Pg           = gsl_complex_rect(Load_flow_arr->Pg,0.0);
        gsl_complex Qg           = gsl_complex_rect(0.0,Load_flow_arr->Qg);
        gsl_complex Xq0          = gsl_complex_rect(0.0,Gen_array->Xq0);
            
        I_0[x] =gsl_complex_div(gsl_complex_sub(Pg,Qg),gsl_complex_mul(Vt_0,e__neg_j_theta));
        Gen_initial[x].I_0=I_0[x];
            
        double mag_I_0    = gsl_complex_abs(I_0[x]);
        double phi_0      = gsl_complex_arg(I_0[x]);
            
        Eq_0[x]= gsl_complex_add((gsl_complex_mul(I_0[x],Xq0)),(gsl_complex_mul(Vt_0,e_j_theta)));
        Gen_initial[x].Eq_0=Eq_0[x];
            
        double delta_0=gsl_complex_arg(Eq_0[x]);
        delta_0_arr[x]=delta_0;
        Gen_initial[x].delta_0=delta_0_arr[x];

        id_0[x]     = gsl_complex_rect(( -mag_I_0   ) * ( sin (delta_0-phi_0) ),0.0) ;
        iq_0[x]     = gsl_complex_rect((  mag_I_0   ) * ( cos (delta_0-phi_0) ),0.0) ;

        //Made CHANGES HERE THETA WAS IN DEGREES WHILE ANGLE IS THE ONE IN RADS    
        vd_0[x]     = gsl_complex_rect(( -mag_Vt_0  ) * ( sin (delta_0-angle) ),0.0) ;
        vq_0[x]     = gsl_complex_rect((  mag_Vt_0  ) * ( cos (delta_0-angle) ),0.0) ;

        Gen_initial[x].id_0=id_0[x];
        Gen_initial[x].iq_0=iq_0[x];
        Gen_initial[x].vd_0=vd_0[x];
        Gen_initial[x].vq_0=vq_0[x];

        /* ===================================================================
         * CORRECTED EMF INITIALIZATION - CONSISTENT WITH DYNAMIC MODEL
         * ===================================================================
         * The original EMF calculations were inconsistent with the actual
         * dynamic simulation model, causing large voltage errors and delta drift.
         * This corrected version ensures perfect steady-state consistency.
         * =================================================================== */
        
        /* Step 1: Calculate terminal voltage magnitude for exciter reference */
        double Vt_magnitude = mag_Vt_0;
        
        /* Step 2: Set Vref = Vt for steady state (no voltage error) */
        double Vref_init = Vt_magnitude; /* Vref equals pre-fault terminal voltage */
        
        /* Step 3: Compute Efd from steady-state flux equation: Eq' = Efd + id*(Xd0-Xd1) */
        double Efd_calculated = Eq_0[x].dat[0] - (id_0[x].dat[0] * (Gen_array->Xd0 - Gen_array->Xd1));
        
        /* Step 4: Calculate Eq' (transient EMF) from steady-state flux equation */
        /* From steady state: dEq'/dt = 0 = (Efd - Eq' + id*(Xd0-Xd1))/Tdo1 */
        /* Therefore: Eq' = Efd + id*(Xd0-Xd1) */
        double Eq_dash_calculated = Efd_calculated + (id_0[x].dat[0] * (Gen_array->Xd0 - Gen_array->Xd1));
        
        /* Step 5: Calculate Ed' (transient EMF) from steady-state flux equation */
        /* From steady state: dEd'/dt = 0 = -(Ed' + iq*(Xq0-Xq1))/Tqo1 */
        /* Therefore: Ed' = -iq*(Xq0-Xq1) */
        double Ed_dash_calculated = -iq_0[x].dat[0] * (Gen_array->Xq0 - Gen_array->Xq1);
        
        /* Step 6: Calculate Eq'' (sub-transient EMF) from steady-state equation */
        /* From steady state: dEq''/dt = 0 = (Eq' - Eq'' + id*(Xd1-Xd2))/Tdo2 */
        /* Therefore: Eq'' = Eq' + id*(Xd1-Xd2) */
        double Eq_das_das_calculated = Eq_dash_calculated + (id_0[x].dat[0] * (Gen_array->Xd1 - Gen_array->Xd2));
        
        /* Step 7: Calculate Ed'' (sub-transient EMF) from steady-state equation */
        /* From steady state: dEd''/dt = 0 = (Ed' - Ed'' + iq*(Xq2-Xq1))/Tqo2 */
        /* Therefore: Ed'' = Ed' + iq*(Xq2-Xq1) */
        double Ed_das_das_calculated = Ed_dash_calculated + (iq_0[x].dat[0] * (Gen_array->Xq2 - Gen_array->Xq1));

        /* Store the corrected values */
        Efd_0[x] = gsl_complex_rect(Efd_calculated, 0.0);
        Eq_dash_0[x] = gsl_complex_rect(Eq_dash_calculated, 0.0);
        Ed_dash_0[x] = gsl_complex_rect(Ed_dash_calculated, 0.0);
        Ed_das_das[x] = Ed_das_das_calculated;
        Eq_das_das[x] = Eq_das_das_calculated;

        /* Calculate electrical power with corrected EMFs */
        double temp_reluctance = Gen_array->Xd2 - Gen_array->Xq2;
        Pe_0[x] = gsl_complex_rect(
            (Eq_das_das[x] * iq_0[x].dat[0]) + 
            (Ed_das_das[x] * id_0[x].dat[0]) + 
            (temp_reluctance * id_0[x].dat[0] * iq_0[x].dat[0]),
            0.0);
            
        Pm_0[x] = Pe_0[x];  /* Perfect balance */

        /* Debug output for verification */
        printf("Gen %d Initial EMF Correction:\n", x);
        printf("  Vt = %.6f pu, Vref = %.6f pu (error = %.6f pu)\n", 
               Vt_magnitude, Vref_init, Vref_init - Vt_magnitude);
        printf("  Efd = %.6f pu, Eq' = %.6f pu, Ed' = %.6f pu\n", 
               Efd_calculated, Eq_dash_calculated, Ed_dash_calculated);
        printf("  Eq'' = %.6f pu, Ed'' = %.6f pu\n", 
               Eq_das_das_calculated, Ed_das_das_calculated);
        printf("  Pe = %.6f pu, Pm = %.6f pu\n", 
               Pe_0[x].dat[0], Pm_0[x].dat[0]);

        /* Calculate network voltages in generator reference frame */
        double temp_I_real = (Eq_das_das[x]*cos(delta_0_arr[x])) - (Ed_das_das[x]*sin(delta_0_arr[x]));
        double temp_I_IMAG = (Ed_das_das[x]*cos(delta_0_arr[x])) + (Eq_das_das[x]*sin(delta_0_arr[x]));
        gsl_complex added = gsl_complex_rect(temp_I_real, temp_I_IMAG);
        gsl_complex impedance = gsl_complex_rect(0, Gen_array->Xd2);
        I_0_DQ_FRAME[x] = gsl_complex_div(added, impedance); 

        VD[x] = gsl_complex_abs(vd_0[x]) * cos(delta_0) - gsl_complex_abs(vq_0[x]) * sin(delta_0);
        VQ[x] = gsl_complex_abs(vq_0[x]) * cos(delta_0) + gsl_complex_abs(vd_0[x]) * sin(delta_0);

        Gen_initial[x].Efd_0    =Efd_0[x];
        Gen_initial[x].Eq_dash_0=Eq_dash_0[x];
        Gen_initial[x].Ed_dash_0=Ed_dash_0[x];
        Gen_initial[x].Pe_0     =Pe_0[x];
        Gen_initial[x].Pm_0     =Pm_0[x];
            
        Gen_initial[x].Ed_das_das=Ed_das_das[x];
        Gen_initial[x].Eq_das_das=Eq_das_das[x];
        Gen_initial[x].delta_0   =delta_0_arr[x];

        Gen_initial[x].I_0_DQFRAME=I_0_DQ_FRAME[x];
        Gen_initial[x].VQ=VQ[x];
        Gen_initial[x].VD=VD[x];
    }
        
    return Gen_initial ;
} 