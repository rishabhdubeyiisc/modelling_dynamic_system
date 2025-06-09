#include    <stdio.h>
#include    <math.h>
#include    <stdlib.h>
#include    <sys/stat.h>
#include    <sys/types.h>

// HEADER FILES OF GSL
#include    <gsl/gsl_complex.h>
#include    <gsl/gsl_complex_math.h> 
// HEADER FILES I CREATED

#include "read_header.h"
#include "data_types.h"
#ifndef NG_MAX
#define NG_MAX 20
#endif

/* forward declaration of the new simplified solver */
void partitioned_solver_sm(
        COMB_STRUCT All_data,
        INITIAL   *Initial_state,
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

// constants

int main()
{   
    //***************************************************//
    
    // double del_t        =(0.0002);  // 0.2 ms is step size
    double del_t        =(0.002);  // 0.2 ms is step size
    double FAULT_CYCLES =800  ; //cycles
    double not_used =0.0; 
    
    double VREF_CHANGE_TIME       = 10; // Changed from 30s to 10s for faster testing
    double MECH_POWER_CHANGE_TIME = 1.5 * 60; // MINUTES -> SECONDS
    double CLEAR_DISTURBANCES     = 2.5 * 60; // MINUTES -> SECONDS
    double FAULT_INITIATION_TIME  = 3.0 * 60; // MINUTES -> SECONDS
    int REMOVE_DISTURBANCES       = 0  ; // TRUE
    int END_TIME_MINUTES          = 10 ; // MINUTES
    int END_TIME                  = END_TIME_MINUTES * 60; // Convert to seconds
    // int DISTURBANCE_CYCLES;
    
    int Fault_Do    =1; //FLAGS
    int Fault_stop  =1; //FLAGS
    //**************************************************//

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
    COMB_STRUCT Net_data_main  =  Reader_fn();
    //creating a variable of suitable dimention blocks to get the initial values as initial values has many things nside it
    INITIAL * Initial_state_main = Gen_data(Net_data_main);
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
    

    printf("\nSimulator will write per-generator CSV files to sim/ directory\n");

    /* Ask user for Vref step parameters */
    double vref_step_time  = 200.0;  /* default 100 s */
    double vref_step_delta = 0.05;   /* default +0.05 pu */
    double pm_step_time    = 200.0;  /* default same as Vref */
    double pm_step_delta   = 0.0;    /* default 0 (no change) */
    int    target_g        = 0;      /* default generator index */

    printf("Enter generator index, Vref Δ (pu), Pm Δ (pu) separated by space (e.g. 1 0.05 0.02): ");
    fflush(stdout);
    if (scanf("%d %lf %lf", &target_g, &vref_step_delta, &pm_step_delta) != 3) {
        printf("Using defaults: gen=0, Vref Δ = 0.05 pu, Pm Δ = 0.0 pu\n");
        target_g = 0; vref_step_delta = 0.05; pm_step_delta = 0.0;
    }

    printf("Enter step time in seconds (default 100): ");
    fflush(stdout);
    if (scanf("%lf", &vref_step_time) != 1) { vref_step_time = 100.0; }

    pm_step_time = vref_step_time; /* same moment for both steps */

    /* Show current reference values */
    printf("Selected Generator %d | Initial Vt ≈ %.4f pu | Initial Pm ≈ %.4f pu\n",
           target_g,
           sqrt(pow(Initial_state_main[target_g].vq_0.dat[0],2)+pow(Initial_state_main[target_g].vd_0.dat[0],2)),
           Initial_state_main[target_g].Pm_0.dat[0]);
    printf("Scheduled: Vref Δ = %+0.4f pu, Pm Δ = %+0.4f pu at t = %.1f s\n",
           vref_step_delta, pm_step_delta, vref_step_time);

    // getting Y_BUS
    
    Y_STRUCT Y_complex_aug_main = Y_BUS_AUG(Net_data_main);
    
    //check the result
  
    double ** Y_Aug_sp_main   = Y_spitter (Y_complex_aug_main,Net_data_main);
    int num_buses = Net_data_main.constants[0].LOAD_FLOW_number;
    double ** Z_AUG           = Y_spitter (Y_complex_aug_main,Net_data_main);
    inv_mat(Y_Aug_sp_main, 2 * num_buses, Z_AUG);

    /* ------------------------------------------------------------------
     *  Fault configuration                                             
     * ------------------------------------------------------------------*/
    double **Z_AUG_fault = NULL;
    int    fault_enabled = 0;
    double fault_start   = 0.0;
    double fault_end     = 0.0;

    printf("\nSimulate fault? (y/n): ");
    fflush(stdout);
    char fault_ans = 'n';
    scanf(" %c", &fault_ans);

    if (fault_ans == 'y' || fault_ans == 'Y') {
        fault_enabled = 1;

        /* show generator-to-bus mapping */
        int Ng_map = Net_data_main.constants[0].NumberofGens;
        int gen_bus[NG_MAX]; /* use same limit as in solver */
        for (int g = 0; g < Ng_map && g < NG_MAX; ++g) {
            gen_bus[g] = Net_data_main.Load_flow_ps[g].bus_number; /* assumption: order matches */
        }

        printf("Available buses (G = generator bus):\n");
        for (int b = 1; b <= num_buses; ++b) {
            int is_gen = -1;
            for (int g = 0; g < Ng_map; ++g)
                if (gen_bus[g] == b) { is_gen = g; break; }
            if (is_gen >= 0)
                printf("  Bus %d  (Gen %d)\n", b, is_gen);
            else
                printf("  Bus %d\n", b);
        }

        printf("\nAvailable lines (annotated with connected generators):\n");
        int num_lines = Net_data_main.constants[0].Number_of_lines;
        for (int i = 0; i < num_lines; ++i) {
            int b1 = Net_data_main.trans_line_para[i].bus1;
            int b2 = Net_data_main.trans_line_para[i].bus2;
            /* mark if either end is a gen bus */
            int g1 = -1, g2 = -1;
            for (int g = 0; g < Ng_map; ++g) {
                if (gen_bus[g] == b1) g1 = g;
                if (gen_bus[g] == b2) g2 = g;
            }
            char tag1[8]="", tag2[8]="";
            if (g1 >= 0) sprintf(tag1," (G%d)", g1);
            if (g2 >= 0) sprintf(tag2," (G%d)", g2);
            printf("  %2d) Bus %d%s — Bus %d%s\n", i + 1, b1, tag1, b2, tag2);
        }

        /* Ask fault type */
        char fault_type = 'l';
        printf("Fault type — terminal bus (t) or line (l): ");
        fflush(stdout);
        scanf(" %c", &fault_type);

        double fault_duration = 0.0;
        printf("Enter fault start time (s) and duration (s): ");
        fflush(stdout);
        scanf("%lf %lf", &fault_start, &fault_duration);
        fault_end = fault_start + fault_duration;

        if (fault_type == 't' || fault_type == 'T') {
            int bus_idx = 1;
            printf("Select bus index (1-%d): ", num_buses);
            fflush(stdout);
            scanf("%d", &bus_idx);

            /* Duplicate healthy admittance and ground the selected bus */
            Y_STRUCT Y_fault = Y_complex_aug_main; /* shallow copy OK */
            /* Ground via huge shunt admittance to force V ≈ 0 */
            Y_fault.MAT[bus_idx - 1][bus_idx - 1].dat[0] = 0.0;
            Y_fault.MAT[bus_idx - 1][bus_idx - 1].dat[1] = 1e6;

            double **Y_AUG_sp_fault = Y_spitter(Y_fault, Net_data_main);
            Z_AUG_fault = Y_spitter(Y_fault, Net_data_main);
            inv_mat(Y_AUG_sp_fault, 2 * num_buses, Z_AUG_fault);

            printf("Configured 3-φ bus fault at Bus %d from %.3f to %.3f s\n", bus_idx, fault_start, fault_end);

        } else { /* default to line fault */
            int line_idx = 1;
            printf("Select line index (1-%d): ", num_lines);
            fflush(stdout);
            scanf("%d", &line_idx);

            Y_STRUCT Y_fault = Y_FAULT_MAKER(line_idx, Net_data_main);
            double **Y_AUG_sp_fault = Y_spitter(Y_fault, Net_data_main);
            Z_AUG_fault = Y_spitter(Y_fault, Net_data_main);
            inv_mat(Y_AUG_sp_fault, 2 * num_buses, Z_AUG_fault);

            int b1 = Net_data_main.trans_line_para[line_idx - 1].bus1;
            int b2 = Net_data_main.trans_line_para[line_idx - 1].bus2;
            printf("Configured line fault %d (Bus %d — Bus %d) from %.3f to %.3f s\n",
                   line_idx, b1, b2, fault_start, fault_end);
        }
    }

    /* use streamlined solver (no fault logic) */
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
    //****************************************************************************//    
    //******************************************************************//
    fclose(fptr_wrte);
    printf("\n\nPartitioned app done\n");
    return 0;    
}



