#include    <stdio.h>
#include    <math.h>
#include    <stdlib.h>

// HEADER FILES OF GSL
#include    <gsl/gsl_complex.h>
#include    <gsl/gsl_complex_math.h> 
// HEADER FILES I CREATED

#include "read_header.h"
#include "data_types.h"


// constants

int main()
{   
    //***************************************************//
    double del_t        =(T_dummy*0.05);
    double FAULT_CYCLES =8.0  ; //cycles
    double not_used =0.0; 
    
    double VREF_CHANGE_TIME       = 0.5; // MINUTES
    double MECH_POWER_CHANGE_TIME = 1.5; // MINUTES
    double CLEAR_DISTURBANCES     = 2.5; // MINUTES
    double FAULT_INITIATION_TIME  = 3.0; // MINUTES
    int REMOVE_DISTURBANCES       = 0  ; // TRUE
    int END_TIME                  = 10 ; // MINUTES
    // int DISTURBANCE_CYCLES;
    
    int Fault_Do    =1; //FLAGS
    int Fault_stop  =1; //FLAGS
    //**************************************************//
    printf("\nWOULD YOU LIKE TO REMOVE THE DISTTURBANCES Y[1]/N[0]\n");
    scanf("%d",&REMOVE_DISTURBANCES);
    //**************************************************//

    FILE *fptr_wrte;

    fptr_wrte = fopen("0_CALCULAT_FILE.txt","w");
    if(fptr_wrte == NULL)
    {
      printf("Error! for write.txt");   
      exit(1);             
    }   
    printf("\nWRITING FILE FOR CALCULATIONS::::::::::::::::::::0_CALCULAT_FILE.txt\n");
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
    

    //**********************************MY_ALGO_MY_ALGO ***************************//
    //**********************************MY_ALGO_MY_ALGO ***************************//
    //**********************************MY_ALGO_MY_ALGO ***************************//
    //**********************************MY_ALGO_MY_ALGO ***************************//
    //**********************************MY_ALGO_MY_ALGO ***************************//

    // FILE * FINAL_TEXT_MY_ALGO;
    // FINAL_TEXT_MY_ALGO = fopen("PLOT_MY_ALGO.csv","w");
    // if ( FINAL_TEXT_MY_ALGO == NULL)
    // {
    //   printf("Error! opening PLOT_MY_ALGO.csv");   
    //   exit(1); 
    // }
    // printf("\nWRITING FILE FOR PLOTTING::::::::::::::::::::::::PLOT_MY_ALGO.csv\n");

    // Y_STRUCT TRIAL_Y_CMPLX_NOT_USEFUL    = Y_MAKER (Net_data_main);
    // Y_STRUCT TRIAL_Y_AUG_MADE            = Y_BUS_AUG_MAKER_NOT_FROM_DATA (TRIAL_Y_CMPLX_NOT_USEFUL,Net_data_main);
    // double ** Y_AUG_SPLIT_TRIAL          = Y_spitter (TRIAL_Y_AUG_MADE,Net_data_main);
    // double ** Z_AUG_SPLIT_TRIAL          = Y_spitter (TRIAL_Y_AUG_MADE,Net_data_main);
    // 
    // inv_mat(Y_AUG_SPLIT_TRIAL,18,Z_AUG_SPLIT_TRIAL);
    // printf("\n Y_BUS_SPLITTED \n");
    // //initializing partitioned approach
    
    // partitioned_solver_new
    //                       (
    //                         Net_data_main,
    //                         Initial_state_main,
    //                         del_t,
    //                         END_TIME,
    //                         not_used,
    //                         FAULT_INITIATION_TIME,
    //                         Fault_Do,
    //                         Fault_stop,
    //                         fptr_wrte,
    //                         Z_AUG_SPLIT_TRIAL,
    //                         Y_AUG_SPLIT_TRIAL,
    //                         TRIAL_Y_AUG_MADE,
    //                         FINAL_TEXT_MY_ALGO,
    //                         VREF_CHANGE_TIME,
    //                         MECH_POWER_CHANGE_TIME,
    //                         CLEAR_DISTURBANCES,
    //                         REMOVE_DISTURBANCES,
    //                         FAULT_CYCLES
    //                       );
    
    // fclose (FINAL_TEXT_MY_ALGO);
    
    //// ******************************SAVED_Y_BUS_SAVED_Y_BUS******************//
    //// ******************************SAVED_Y_BUS_SAVED_Y_BUS******************//
    //// ******************************SAVED_Y_BUS_SAVED_Y_BUS******************//
    //// ******************************SAVED_Y_BUS_SAVED_Y_BUS******************//
    //// ******************************SAVED_Y_BUS_SAVED_Y_BUS******************//

    FILE * FINAL_TEXT_SAVED_BUS;
    FINAL_TEXT_SAVED_BUS = fopen("PLOT_SAVED_BUS.csv","w");
    if ( FINAL_TEXT_SAVED_BUS == NULL)
    {
      printf("Error! opening PLOT_SAVED_BUS.csv");   
      exit(1); 
    }
    printf("\nWRITING FILE FOR PLOTTING::::::::::::::::::::::::PLOT_SAVED_BUS.csv\n");

    // getting Y_BUS
    
    Y_STRUCT Y_complex_aug_main = Y_BUS_AUG(Net_data_main);
    
    //check the result
  
    double ** Y_Aug_sp_main   = Y_spitter (Y_complex_aug_main,Net_data_main);
    double ** Z_AUG           = Y_spitter (Y_complex_aug_main,Net_data_main);
    inv_mat(Y_Aug_sp_main,18,Z_AUG);

    partitioned_solver_new
                          (
                            Net_data_main,
                            Initial_state_main,
                            del_t,
                            END_TIME,
                            not_used,
                            FAULT_INITIATION_TIME,
                            Fault_Do,
                            Fault_stop,
                            fptr_wrte,
                            Z_AUG,
                            Y_Aug_sp_main,
                            Y_complex_aug_main,
                            FINAL_TEXT_SAVED_BUS,
                            VREF_CHANGE_TIME,
                            MECH_POWER_CHANGE_TIME,
                            CLEAR_DISTURBANCES,
                            REMOVE_DISTURBANCES,
                            FAULT_CYCLES
                          );
    //****************************************************************************//    
    fclose (FINAL_TEXT_SAVED_BUS);
    //******************************************************************//
    fclose(fptr_wrte);
    printf("\n\nPartitioned app done\n");
    return 0;    
}



