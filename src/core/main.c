#include    <stdio.h>
#include    <math.h>
#include    <stdlib.h>
#include    <sys/stat.h>
#include    <sys/types.h>

// HEADER FILES OF GSL
#include    <gsl/gsl_complex.h>
#include    <gsl/gsl_complex_math.h> 
// HEADER FILES I CREATED

#include "../../include/api/power_system.h"
#include "../../include/core/simulator.h"

// constants

int main()
{   
    //***************************************************//
    
    // double del_t        =(0.0002);  // 0.2 ms is step size works very good 
    
    double del_t = 0.0002;  // 2 ms time step
    int END_TIME_MINUTES = 10;     // MINUTES  
    int END_TIME = END_TIME_MINUTES * 60; // Convert to seconds

    FILE *fptr_wrte;

    /* Create simulation output directory */
    mkdir("sim", 0755);

    fptr_wrte = fopen("sim/sim_prints.txt","w");
    if(fptr_wrte == NULL)
    {
      printf("Error! opening sim/sim_prints.txt");   
      exit(1);             
    }   
    printf("\nWRITING FILE FOR CALCULATIONS::::::::::::::::::::sim/sim_prints.txt\n");
    char text_head[]= " TIME      DELTA     SLIP      Efd      Ed2       Eq2     Ed1      Eq1     [0].VQ   [0].VD   vq[0]     vd[0]    Vt[0]     Deg     Error[0]   iq       id       Elec";
    fprintf(fptr_wrte, " %s \n ",text_head);

    //this is to read data from the file 
    NetworkData Net_data_main  =  loadSystemData();
    //creating a variable of suitable dimention blocks to get the initial values as initial values has many things nside it
    InitialConditions * Initial_state_main = initializeGenerator(Net_data_main);
    printf("\n INITIAL_STAGE_DONE \n");
    fprintf(fptr_wrte," \n INITIAL VALUES \n ");
    fprintf(fptr_wrte," %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  \n ",
                      0.0,
                      Initial_state_main[0].delta_0 ,
                      0.0,
                      Initial_state_main[0].Efd_0.dat[0] ,
                      Initial_state_main[0].Ed_das_das ,
                      Initial_state_main[0].Eq_das_das,
                      Initial_state_main[0].Ed_dash_0.dat[0],
                      Initial_state_main[0].Eq_dash_0.dat[0],
                      Initial_state_main[0].VQ,
                      Initial_state_main[0].VD,
                      Initial_state_main[0].vq_0.dat[0],
                      Initial_state_main[0].vd_0.dat[0],
                      sqrt(pow(Initial_state_main[0].vq_0.dat[0],2)+pow(Initial_state_main[0].vd_0.dat[0],2)),
                      Initial_state_main[0].delta_0*(180/PI),
                      0.0,
                      Initial_state_main[0].iq_0.dat[0],
                      Initial_state_main[0].id_0.dat[0] 
                      );  
    
    fclose(fptr_wrte);

    printf("\nSimulator will write per-generator CSV files to sim/ directory\n");

    /* Get simulation parameters from user */
    double vref_step_time, vref_step_delta, pm_step_time, pm_step_delta;
    int target_g;
    
    get_simulation_parameters(Net_data_main.constants[0].NumberofGens,
                             &target_g, &vref_step_delta, &pm_step_delta,
                             &vref_step_time, &pm_step_time, Initial_state_main);

    // getting Y_BUS
    
    AdmittanceMatrix Y_complex_aug_main = createAugmentedAdmittanceMatrix(Net_data_main);
    
    //check the result
  
    double ** Y_Aug_sp_main   = convertComplexToReal (Y_complex_aug_main,Net_data_main);
    int num_buses = Net_data_main.constants[0].LOAD_FLOW_number;
    double ** Z_AUG           = convertComplexToReal (Y_complex_aug_main,Net_data_main);
    invertMatrix(Y_Aug_sp_main, 2 * num_buses, Z_AUG);

    /* Configure fault simulation */
    double **Z_AUG_fault = NULL;
    double fault_start = 0.0;
    double fault_end = 0.0;
    int fault_enabled = configure_fault_simulation(Net_data_main, Y_complex_aug_main, 
                                                   &Z_AUG_fault, &fault_start, &fault_end);

    /* Debug matrix prints removed for cleaner output */

    /* Choose solver method */
    printf("\nSelect solver method:\n");
    printf("  1) Partitioned solver (fast, production) \n");
    printf("  2) Jacobian solver   (simultaneous DAE)\n");
    printf("Choice (1 or 2) [default 1]: ");
    fflush(stdout);
    
    int solver_choice = 1;   /* default → partitioned */
    if (scanf("%d", &solver_choice) != 1) {
        solver_choice = 1;   /* user just hit <Enter> */
    }
    /* flush the rest of the line to avoid dangling newline before next scanf/fgets */
    int c; while ((c = getchar()) != '\n' && c != EOF);
    
    if (solver_choice != 2) solver_choice = 1;  /* guard */
    
    if (solver_choice == 1) {
        printf("✅ Running partitioned solver...\n");
        partitioned_solver_sm(
            Net_data_main,
            Initial_state_main,
            del_t,
            END_TIME,
            Z_AUG,
            Z_AUG_fault,
            fault_enabled,
            fault_start,
            fault_end,
            target_g,
            vref_step_time,
            vref_step_delta,
            pm_step_time,
            pm_step_delta
        );
      printf("\n\nPartitioned app done\n");
    } else {
        printf("✅ Running Jacobian solver...\n");
        jacobian_solver_sm(
            Net_data_main,
            Initial_state_main,
            del_t,
            END_TIME,
            Z_AUG,
            Z_AUG_fault,
            fault_enabled,
            fault_start,
            fault_end,
            target_g,
            vref_step_time,
            vref_step_delta,
            pm_step_time,
            pm_step_delta
        );
      printf("\n\nJacobian app done\n");
    }
    //****************************************************************************//    
    //******************************************************************//
    return 0;    
}



