#include    <stdio.h>
#include    <math.h>
#include    <stdlib.h>

// HEADER FILES OF GSL
#include    <gsl/gsl_complex.h>
#include    <gsl/gsl_complex_math.h>

// HEADER FILES I CREATED

#include    "read_header.h"
#include    "data_types.h"

double Power_elc(double Ed2,double id,double Eq2,double iq,double Xd2,double Xq2)
{       
        double power_elec;
        double term1 = Ed2*id;
        double term2 = Eq2*iq;
        double term3 = (Xd2-Xq2)*(id)*(iq);
        
        power_elec      =term1+term2+term3;
        return power_elec;
}
//*******************************************//
void partitioned_solver_new
        (
        COMB_STRUCT All_data,
        INITIAL * Initial_state,
        double del_t,
        double END_TIME,
        double not_used,
        double FAULT_INITIATION_TIME,
        int FAULT_DO_FLAG,
        int FAULT_STOP_FLAG,
        FILE *fptr,
        double ** Z_AUG_healthy,
        double ** Y_AUG_healthy,
        Y_STRUCT Y_AUG_complex_healthy,
        FILE * FINAL_FILE,
        double VREF_CHANGE_TIME,
        double MECH_POWER_ALTERATION_TIME,
        double DISTURBANCES_CLEARING_TIME,
        int REMOVE_DISTURBANCES_FLAG,
        double FAULTED_CYCLES
        )
{   
    double PERIOD = (pow(60,-1));
    double FAULTED_TIME = FAULTED_CYCLES*(PERIOD);

    int Number_of_generators=All_data.constants[0].NumberofGens;
    int Number_of_buses     =All_data.constants[0].LOAD_FLOW_number;

    double D_by_m[Number_of_generators];
    D_by_m[0]=0.1;
    D_by_m[1]=0.2;  
    D_by_m[2]=0.3;
    
    double Xd2_local[Number_of_generators];
    double Xq2_local[Number_of_generators];
        
    for (int i = 0; i < Number_of_generators; i++)
    {
        Xd2_local[i]=All_data.Generator_ps[i].Xd2;
        Xq2_local[i]=All_data.Generator_ps[i].Xq2;
    }

    double Vref[Number_of_generators];
    double Vt[Number_of_generators];
    double * ptr_Vt;
    ptr_Vt=&Vt[0];

    double iq[Number_of_generators];
    double * ptr_iq;
    ptr_iq = &iq[0];
    
    double id[Number_of_generators];
    double * ptr_id;
    ptr_id = &id[0];

    double vd[Number_of_generators];
    double vq[Number_of_generators];

    double mech_power[Number_of_generators];
    double Elec_power[Number_of_generators];

    STATES X_elr[Number_of_generators];
    STATES * pointer_to_X_elr;
    pointer_to_X_elr=&X_elr[0];

    STATES X_vector[Number_of_generators];
    STATES * pointer_X_VECTOR;
    pointer_X_VECTOR = &X_vector[0];    
        
    DIFF_STATES_VECTOR F_vector_old[Number_of_generators];
    DIFF_STATES_VECTOR  * pointer_F_VECTOR_old ;
    pointer_F_VECTOR_old = &F_vector_old[0];

    DIFF_STATES_VECTOR F_vector_new[Number_of_generators];
    DIFF_STATES_VECTOR  * pointer_F_VECTOR_new ;
    pointer_F_VECTOR_new = &F_vector_new[0];

    DIFF_STATES_VECTOR F_Sum[Number_of_generators];
    DIFF_STATES_VECTOR * pointer_to_Fsum;
    pointer_to_Fsum = &F_Sum[0];

    NW_STATES  Nw_vector[Number_of_buses];    

    //initialize this vector for the other remaing buses
        
    NW_STATES * pointer_Nw_vector;
    pointer_Nw_vector=&Nw_vector[0];

    double proportional_error[Number_of_generators];
    //8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

    int initial_status=0;

    int vref_change_Status=0;
    int asked_vref =0; 
    double Vref_user;
    int Exiter_change_location ;

    int Prime_mover_change_status = 0 ;
    int Prime_mover_change_location;
    int asked_Prime_mover=0;
    double Prime_mover_user_input;
    
    int FAULT_ANSWER=99; 
    int FAULT_STATUS=0;
    int FAULT_LOCATION_BUS_NUMBER;
    int FAULT_ASKED=0;
    char FAULT_TYPE = 'F';
    int FAULTED_LINE_SET;
    char TERMINAL = 'T';
    char LINES    = 'L';
    Y_STRUCT Y_faulted;
    int Augment_calculated=0;
    double ** Y_AUG_fault_split;
    double ** Z_AUG_fault_split;

    //8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888    
    fprintf(FINAL_FILE,"time,delta,slip,Efd,Ed2,Eq2,Ed1,Eq1,VQ_0,VD_0,Vref,Vt,EROR,iq_0,id_0,mech_power,elec_power\n"); 
    
    //start loop 
    for (double t = 0; t < END_TIME; t=t+del_t)
    {   //ptr_to_All_data=&All_data;
        if (t==0 && initial_status == 0)   // Fold=Fnew and start i.e Xold in Fnew put and calc or just take 2 times Fold
        {
                for (int i = 0; i < Number_of_generators ; i++)//setting the initial values of v ref
                {  
                
                        Vt[i]= sqrt(    pow(gsl_complex_abs(Initial_state[i].vd_0),2)+
                                pow(gsl_complex_abs(Initial_state[i].vq_0),2)
                        );
                
                
                        Vref[i]=(Vt[i])+((1/All_data.exciter_ps[i].ka) * (Initial_state[i].Efd_0.dat[0]));

                        proportional_error[i]=(Vref[i]-Vt[i]);

                        vd[i]=(Initial_state[i].vd_0.dat[0]);
                        vq[i]=(Initial_state[i].vq_0.dat[0]);

                        X_vector[i].Ed_das=(Initial_state[i].Ed_dash_0.dat[0]);
                        X_vector[i].Eq_das=(Initial_state[i].Eq_dash_0.dat[0]);
                        X_vector[i].Ed_das_das=Initial_state[i].Ed_das_das;
                        X_vector[i].Eq_das_das=Initial_state[i].Eq_das_das;
                        X_vector[i].delta=Initial_state[i].delta_0;
                        X_vector[i].slip=0;
                        X_vector[i].Efd=(Initial_state[i].Efd_0.dat[0]);
                        X_vector[i].E_dummy=(Xd2_local[i]-Xq2_local[i])*(Initial_state[i].iq_0.dat[0]); 

                        id[i]=(Initial_state[i].id_0.dat[0]) ;
                        iq[i]=(Initial_state[i].iq_0.dat[0]) ;
                        mech_power[i]=(Initial_state[i].Pm_0.dat[0]) ;
                        Elec_power[i]=(Initial_state[i].Pe_0.dat[0]) ;
                        Nw_vector[i].VQ=Initial_state[i].VQ;
                        Nw_vector[i].VD=Initial_state[i].VD;
                }
                //******************************INITIAL_VALUES_TAKEN*******************************************************************//
                //***********************NETWORK SOLVED FROM INITIAL VALUES *********************************************//
                Network_solver( All_data,
                                Number_of_generators,
                                Number_of_buses,
                                pointer_Nw_vector,
                                pointer_X_VECTOR,
                                Z_AUG_healthy,
                                ptr_iq,
                                ptr_id,
                                ptr_Vt);
                //******************************
                pointer_Nw_vector=&Nw_vector[0];
                pointer_X_VECTOR=&X_vector[0];
                ptr_iq = & iq[0];
                ptr_id = & id[0];
                //**************************************NETWRK SOLVED ***********************************************************************************// 
                for (int i = 0; i < Number_of_generators; i++)
                {
                        Vt[i]= sqrt (pow(Nw_vector[i].VQ,2)+pow(Nw_vector[i].VD,2));
                        vq[i]=Nw_vector[i].VQ*cos(X_vector[i].delta)+Nw_vector[i].VD*sin(X_vector[i].delta);
                        vd[i]=Nw_vector[i].VD*cos(X_vector[i].delta)-Nw_vector[i].VQ*sin(X_vector[i].delta);
                } 
                //***************************Updated Vt from new netwrok states*****************************************//
                //*****************************SOLVE FOR EFD ***************************************************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_X_VECTOR++)
                {       
                        proportional_error[i]=Vref[i]-Vt[i];
                        Efd_E_dummy_solver(   i,
                                        All_data,
                                        pointer_X_VECTOR,
                                        proportional_error[i],
                                        del_t,
                                        iq[i],
                                        Xq2_local[i],
                                        Xd2_local[i]);
                }
                pointer_X_VECTOR = & X_vector[0];
                // **************************CALCULATED EFD HERE ****************************************//
                //****************************CURRENT_UPDATED******************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        ptr_iq++,
                                                        ptr_id++)
                {
                        
                        current_update_new(X_vector[i],
                                        vd[i],
                                        vq[i],
                                        ptr_id,
                                        ptr_iq,
                                        Xd2_local[i]);
                        /* Keep E_dummy algebraically consistent with the
                           freshly updated iq */
                        X_vector[i].E_dummy = -(Xq2_local[i] - Xd2_local[i]) * iq[i];
                }
                ptr_iq=&iq[0];
                ptr_id=&id[0];
                
                //****************************CURRENT_UPDATED******************************************//
                //****************************CURRENT_UPDATED******************************************//
                
                //****************************UPDATING_POWER*****************************************//
                for (int i = 0; i < Number_of_generators; i++)
                {
                        Elec_power[i]=Power_elc(X_vector[i].Ed_das_das,
                                                id[i],
                                                X_vector[i].Eq_das_das,
                                                iq[i],
                                                Xd2_local[i],
                                                Xq2_local[i]);
                }

                //****************************UPDATING_POWER*****************************************//
                //****************************UPDATING_POWER*****************************************//
      
                for (int i = 0; i < Number_of_generators;i++,
                                                         pointer_F_VECTOR_old++)
                {
                        F_VECTOR_CALC_single(i,
                                        pointer_F_VECTOR_old,
                                        All_data,
                                        X_vector[i],
                                        D_by_m[i],
                                        id[i],
                                        iq[i],
                                        mech_power[i],
                                        Elec_power[i]);
                }
                pointer_F_VECTOR_old = & F_vector_old[0];
                //**************************************************************************************//
                //**************************************************************************************//
                for (int i = 0; i < Number_of_generators; ++i,
                                                        pointer_X_VECTOR++,
                                                        pointer_F_VECTOR_old++)
                {
                        /* 1. start from current state */
                        X_elr[i] = *pointer_X_VECTOR;

                        /* 2. predictor */
                        Euler_forward_wo_Exiter(
                                                 &X_elr[i],
                                                 pointer_X_VECTOR,
                                                 pointer_F_VECTOR_old,
                                                 del_t
                                                 );
                }
                pointer_to_X_elr=&X_elr[0];                
                pointer_X_VECTOR = & X_vector[0];
                pointer_F_VECTOR_old= &F_vector_old[0];

                for (int i = 0; i < Number_of_generators; 
                                                        i++,
                                                        pointer_F_VECTOR_new++)
                {
                        F_VECTOR_CALC_single(i,
                                        pointer_F_VECTOR_new,
                                        All_data,
                                        X_elr[i],
                                        D_by_m[i],
                                        id[i],
                                        iq[i],
                                        mech_power[i],
                                        Elec_power[i]);
                }
                pointer_F_VECTOR_new = & F_vector_new[0];
                //f_new_calculated 
                //888888888888888888888888888888888888888
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_to_Fsum++,
                                                        pointer_F_VECTOR_old++,
                                                        pointer_F_VECTOR_new++,
                                                        pointer_X_VECTOR++)
                {
                        Summer_diffential(      pointer_to_Fsum,
                                                pointer_F_VECTOR_old,
                                                pointer_F_VECTOR_new);
                        Sum_Xold_F_rets_Xnew
                                        (       pointer_X_VECTOR,
                                                pointer_to_Fsum,
                                                del_t);
                }
                pointer_to_Fsum         =       & F_Sum[0];
                pointer_F_VECTOR_old    =       &F_vector_old[0];
                pointer_F_VECTOR_new    =       &F_vector_new[0];
                pointer_X_VECTOR        =       &X_vector[0];

                //************** DONE TRAP 1st time *****************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_F_VECTOR_old++,
                                                        pointer_F_VECTOR_new++)
                {
                        copy_mak_set_ptr_zer(pointer_F_VECTOR_new,pointer_F_VECTOR_old);
                }
                pointer_F_VECTOR_old=&F_vector_old[0];
                pointer_F_VECTOR_new=&F_vector_new[0];
      
                initial_status = 1;
        }
              
        //*************************************DOING DISTURBANCES *************************//
        if ( t < DISTURBANCES_CLEARING_TIME && REMOVE_DISTURBANCES_FLAG == 0 )
        {      
                if (REMOVE_DISTURBANCES_FLAG == 0 && initial_status ==1)
                {
                        printf("\nDISTURBANCES ARE NOT GOING TO BE REMOVED\n ");
                        printf("*******************************************\n");
                        printf("\nDISTURBANCES ARE NOT GOING TO BE REMOVED\n ");
                        printf("*******************************************\n");
                        printf("\nDISTURBANCES ARE NOT GOING TO BE REMOVED\n ");
                        printf("*******************************************\n");
                        printf("\nDISTURBANCES ARE NOT GOING TO BE REMOVED\n ");

                        initial_status =2;
                }
                
                if (VREF_CHANGE_TIME<t && asked_vref == 0  )
                {       
                        int answer_vref=99; 
                        printf("\n\nYOU WANT TO CHANGE Vref Y[1],N[0]\n");
                        scanf("%d",&answer_vref);
                        
                        if (answer_vref == 1 )
                        {       
                                printf("LOCATION OF WHICH GENERATOR YOU WANT TO CHANGE 0 1 2\n");
                                scanf("%d",&Exiter_change_location);

                                printf("PRESENT VALUE OF VREF AT LOCATION is ::: %lf \n",Vref[Exiter_change_location]);
                                printf("GIVE VALUE OF VREF\n");
                                scanf("%lf",&Vref_user);
                                
                                vref_change_Status = 1;

                        }
                        if (answer_vref == 0)
                        {
                               printf("NOT CHANGING VREF\n"); 
                        }
                        asked_vref=1;                        
                               
                }

                if (MECH_POWER_ALTERATION_TIME<t && asked_Prime_mover == 0 )
                {
                        int answer_power=99; 
                        printf("\n\nYOU WANT TO CHANGE PRIME MOVER INPUT Y[1],N[0]\n");
                        scanf("%d",&answer_power);
                        
                        if (answer_power == 1 )
                        {       
                                printf("LOCATION OF GENERATOR YOU WANT TO CHANGE PRIME MOVER POWER :::: 0, 1, 2\n");
                                scanf("%d",&Prime_mover_change_location);

                                printf("PRESENT VALUE OF MECH_POWER AT LOCATION is :::::::::: %lf \n\n",mech_power[Prime_mover_change_location]);
                                printf("GIVE VALUE OF MECH_POWER\n");
                                scanf("%lf",&Prime_mover_user_input);
                                
                                Prime_mover_change_status = 1;

                        }
                        if (answer_power == 0)
                        {
                               printf("NOT CHANGING PRIME MOVER POWER \n"); 
                        }
                        asked_Prime_mover=1;              
                }
                if (vref_change_Status == 1)
                {
                        Vref[Exiter_change_location]=Vref_user;
                }
                if (Prime_mover_change_status == 1)
                {
                        mech_power[Prime_mover_change_location]=Prime_mover_user_input;
                }

                
                //***********************PASTE HERE ***********************************************//
                //***********************PASTE HERE ***********************************************//
                //***********************PASTE HERE ***********************************************//
                //***********************PASTE HERE ***********************************************//
                //***********************PASTE HERE ***********************************************//
                //***********************NETWORK SOLVED FROM INITIAL VALUES *********************************************//
                Network_solver( All_data,
                                Number_of_generators,
                                Number_of_buses,
                                pointer_Nw_vector,
                                pointer_X_VECTOR,
                                Z_AUG_healthy,
                                ptr_iq,
                                ptr_id,
                                ptr_Vt);
                //******************************
                pointer_Nw_vector=&Nw_vector[0];
                pointer_X_VECTOR=&X_vector[0];
                ptr_iq = & iq[0];
                ptr_id = & id[0];

                //**************************************NETWRK SOLVED ***********************************************************************************// 
                for (int i = 0; i < Number_of_generators; i++)
                {
                        Vt[i]= sqrt (pow(Nw_vector[i].VQ,2)+pow(Nw_vector[i].VD,2));
                        vq[i]=Nw_vector[i].VQ*cos(X_vector[i].delta)+Nw_vector[i].VD*sin(X_vector[i].delta);
                        vd[i]=Nw_vector[i].VD*cos(X_vector[i].delta)-Nw_vector[i].VQ*sin(X_vector[i].delta);
                } 
                //***************************Updated Vt from new netwrok states*****************************************//
                //*****************************SOLVE FOR EFD ***************************************************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_X_VECTOR++)
                {       
                        proportional_error[i]=Vref[i]-Vt[i];
                        Efd_E_dummy_solver(   i,
                                        All_data,
                                        pointer_X_VECTOR,
                                        proportional_error[i],
                                        del_t,
                                        iq[i],
                                        Xq2_local[i],
                                        Xd2_local[i]);
                }
                pointer_X_VECTOR = & X_vector[0];
                // **************************CALCULATED EFD HERE ****************************************//
                //****************************CURRENT_UPDATED******************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        ptr_iq++,
                                                        ptr_id++)
                {
                        
                        current_update_new(X_vector[i],
                                        vd[i],
                                        vq[i],
                                        ptr_id,
                                        ptr_iq,
                                        Xd2_local[i]);
                        /* Keep E_dummy algebraically consistent with the
                           freshly updated iq */
                        X_vector[i].E_dummy = -(Xq2_local[i] - Xd2_local[i]) * iq[i];
                }
                ptr_iq=&iq[0];
                ptr_id=&id[0];
                
                //****************************CURRENT_UPDATED******************************************//
                //****************************CURRENT_UPDATED******************************************//
                
                //****************************UPDATING_POWER*****************************************//
                for (int i = 0; i < Number_of_generators; i++)
                {
                        Elec_power[i]=Power_elc(X_vector[i].Ed_das_das,
                                                id[i],
                                                X_vector[i].Eq_das_das,
                                                iq[i],
                                                Xd2_local[i],
                                                Xq2_local[i]);
                }


                //****************************UPDATING_POWER*****************************************//
                //****************************UPDATING_POWER*****************************************//
                
                for (int i = 0; i < Number_of_generators;i++,
                                                         pointer_F_VECTOR_old++)
                {
                        F_VECTOR_CALC_single(i,
                                        pointer_F_VECTOR_old,
                                        All_data,
                                        X_vector[i],
                                        D_by_m[i],
                                        id[i],
                                        iq[i],
                                        mech_power[i],
                                        Elec_power[i]);
                }
                pointer_F_VECTOR_old = & F_vector_old[0];
                //**************************************************************************************//
                //**************************************************************************************//
                for (int i = 0; i < Number_of_generators; ++i,
                                                        pointer_X_VECTOR++,
                                                        pointer_F_VECTOR_old++)
                {
                        /* 1. start from current state */
                        X_elr[i] = *pointer_X_VECTOR;

                        /* 2. predictor */
                        Euler_forward_wo_Exiter(&X_elr[i],
                                                 pointer_X_VECTOR,
                                                 pointer_F_VECTOR_old,
                                                 del_t);

                        /* 3. atomic copy-back for this generator */
                        /*
                         * DO NOT copy the Euler-predicted state back to X_vector here.
                         * Keeping X_vector unchanged until after the trapezoidal
                         * corrector preserves the two-step Heun integration scheme
                         * and prevents numerical drift/instability.
                         */
                        /* *pointer_X_VECTOR = X_elr[i]; */
                }     
                pointer_to_X_elr=&X_elr[0];                
                pointer_X_VECTOR = & X_vector[0];
                pointer_F_VECTOR_old= &F_vector_old[0];

                for (int i = 0; i < Number_of_generators; 
                                        ++i,
                                        pointer_F_VECTOR_new++)
                {
                        F_VECTOR_CALC_single(i,
                                        pointer_F_VECTOR_new,
                                        All_data,
                                        X_elr[i],
                                        D_by_m[i],
                                        id[i],
                                        iq[i],
                                        mech_power[i],
                                        Elec_power[i]);
                } 
                pointer_F_VECTOR_new = & F_vector_new[0];
                //f_new_calculated 
                //888888888888888888888888888888888888888
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_to_Fsum++,
                                                        pointer_F_VECTOR_old++,
                                                        pointer_F_VECTOR_new++,
                                                        pointer_X_VECTOR++)
                {
                        Summer_diffential(      pointer_to_Fsum,
                                                pointer_F_VECTOR_old,
                                                pointer_F_VECTOR_new);
                        Sum_Xold_F_rets_Xnew
                                        (       pointer_X_VECTOR,
                                                pointer_to_Fsum,
                                                del_t);
                }
                pointer_to_Fsum         =       & F_Sum[0];
                pointer_F_VECTOR_old    =       &F_vector_old[0];
                pointer_F_VECTOR_new    =       &F_vector_new[0];
                pointer_X_VECTOR        =       &X_vector[0];                
 
                //************** DONE TRAP 1st time *****************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_F_VECTOR_old++,
                                                        pointer_F_VECTOR_new++)
                {
                                copy_mak_set_ptr_zer(pointer_F_VECTOR_new,pointer_F_VECTOR_old);
                }
                pointer_F_VECTOR_old=&F_vector_old[0];
                pointer_F_VECTOR_new=&F_vector_new[0];
       
              
                
        }
        //*************************************FAULT SIMULATION *************************//
        //*************************************FAULT SIMULATION *************************//
        //*************************************FAULT SIMULATION *************************//
        //*************************************FAULT SIMULATION *************************//
        //*************************************FAULT SIMULATION *************************//
        //*************************************FAULT SIMULATION *************************//
        if ( t > 2*DISTURBANCES_CLEARING_TIME && t < (2*DISTURBANCES_CLEARING_TIME + FAULTED_TIME) && REMOVE_DISTURBANCES_FLAG == 1 )
        {                      
                if ( FAULT_ASKED == 0  )
                {       
                       
                        printf("\n\nYOU WANT TO SIMULATE FAULT Y[1],N[0]\n");
                        scanf("%d",&FAULT_ANSWER);
                        
                        if (FAULT_ANSWER == 1 )
                        {       
                                printf("\nWHICH TYPE OF FAULT [T]AT TERMINAL [L]BTW LINES\n ");
                                scanf(" %c", &FAULT_TYPE);

                                if (FAULT_TYPE == TERMINAL )
                                {
                                        printf("BUS number WANT TO SIMULATE FAULT 0 1 2\n");
                                        scanf("%d",&FAULT_LOCATION_BUS_NUMBER);
                                        printf("\nRECORDED BUS NUMBER\n");
                                }
                                if (FAULT_TYPE == LINES)
                                {
                                        printf("LINES NUMBER HAS TO BE GIVEN \n");
                                        printf(" (1 4),::1 \n ") ;                            
                                        printf(" (2 7),::2 \n ") ;                              
                                        printf(" (3 9),::3 \n ") ;                            
                                        printf(" (4 5),::4 \n ") ;                             
                                        printf(" (4 6),::5 \n ") ;                             
                                        printf(" (5 7),::6 \n ") ;                                
                                        printf(" (6 9),::7 \n ") ;                                
                                        printf(" (7 8),::8 \n ") ;                                
                                        printf(" (8 9),::9 \n ") ;
                                        printf("\n\n");                                
                                        scanf(" %d", &FAULTED_LINE_SET);
                                        printf("\nRECORDED BUS NUMBER\n");
                                }
                                FAULT_STATUS = 1;

                        }
                        if (FAULT_ANSWER == 0)
                        {
                               printf("\n NOT DOING FAULT \n"); 
                        }
                        FAULT_ASKED=1;                        
                               
                }
                if (FAULT_STATUS == 1 && Augment_calculated ==0)
                {    
                     //EDIT Y BUS INVERT IT AND APPLY IT
                     if (FAULT_TYPE == TERMINAL)
                     {
                        int Bus_num = FAULT_LOCATION_BUS_NUMBER;
                        Y_AUG_complex_healthy.MAT[Bus_num][Bus_num].dat[0]=0;
                        Y_AUG_complex_healthy.MAT[Bus_num][Bus_num].dat[1]=100;
                        Y_AUG_fault_split   = Y_spitter (Y_AUG_complex_healthy,All_data);
                        Z_AUG_fault_split   = Y_spitter (Y_AUG_complex_healthy,All_data);
                        inv_mat(Y_AUG_fault_split,18,Z_AUG_fault_split);
                     }
                     if (FAULT_TYPE == LINES )
                     {
                        int Line = FAULTED_LINE_SET;
                        Y_faulted = Y_FAULT_MAKER(Line,All_data);
                        Y_AUG_fault_split   = Y_spitter (Y_faulted,All_data);
                        Z_AUG_fault_split   = Y_spitter (Y_faulted,All_data);
                        inv_mat(Y_AUG_fault_split,2*Number_of_buses,Z_AUG_fault_split);
                     }
                     
                Augment_calculated=1;     
                }
                
                //***********************PASTE HERE ***********************************************//
                //***********************PASTE HERE ***********************************************//
                //***********************PASTE HERE ***********************************************//
                //***********************PASTE HERE ***********************************************//
                //***********************PASTE HERE ***********************************************//
                if (FAULT_ANSWER == 0)
                {
                Network_solver( All_data,
                                Number_of_generators,
                                Number_of_buses,
                                pointer_Nw_vector,
                                pointer_X_VECTOR,
                                Z_AUG_healthy,
                                ptr_iq,
                                ptr_id,
                                ptr_Vt);                         
                }
                if (FAULT_ANSWER == 1)
                {
                Network_solver( All_data,
                                Number_of_generators,
                                Number_of_buses,
                                pointer_Nw_vector,
                                pointer_X_VECTOR,
                                Z_AUG_fault_split,
                                ptr_iq,
                                ptr_id,
                                ptr_Vt);
                }
                
                //******************************
                pointer_Nw_vector=&Nw_vector[0];
                pointer_X_VECTOR=&X_vector[0];
                ptr_iq = & iq[0];
                ptr_id = & id[0];
                
                //**************************************NETWRK SOLVED ***********************************************************************************// 
                // 🔧 CRITICAL FIX: Update currents based on predicted states X_elr
                // for correct Heun's method implementation
                double id_elr[Number_of_generators];
                double iq_elr[Number_of_generators];
                
                for (int i = 0; i < Number_of_generators; i++)
                {
                        current_update_new(X_elr[i],
                                        vd[i],
                                        vq[i],
                                        &id_elr[i],
                                        &iq_elr[i],
                                        Xd2_local[i]);
                }
for (int i = 0; i < Number_of_generators; i++)
                {
                        Vt[i]= sqrt (pow(Nw_vector[i].VQ,2)+pow(Nw_vector[i].VD,2));
                        vq[i]=Nw_vector[i].VQ*cos(X_vector[i].delta)+Nw_vector[i].VD*sin(X_vector[i].delta);
                        vd[i]=Nw_vector[i].VD*cos(X_vector[i].delta)-Nw_vector[i].VQ*sin(X_vector[i].delta);
                } 
                //***************************Updated Vt from new netwrok states*****************************************//
                //*****************************SOLVE FOR EFD ***************************************************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_X_VECTOR++)
                {       
                        proportional_error[i]=Vref[i]-Vt[i];
                        Efd_E_dummy_solver(   i,
                                        All_data,
                                        pointer_X_VECTOR,
                                        proportional_error[i],
                                        del_t,
                                        iq[i],
                                        Xq2_local[i],
                                        Xd2_local[i]);
                }
                pointer_X_VECTOR = & X_vector[0];
                // **************************CALCULATED EFD HERE ****************************************//
                //****************************CURRENT_UPDATED******************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        ptr_iq++,
                                                        ptr_id++)
                {
                        
                        current_update_new(X_vector[i],
                                        vd[i],
                                        vq[i],
                                        ptr_id,
                                        ptr_iq,
                                        Xd2_local[i]);
                }
                ptr_iq=&iq[0];
                ptr_id=&id[0];

                //****************************CURRENT_UPDATED******************************************//
                //****************************CURRENT_UPDATED******************************************//
                
                //****************************UPDATING_POWER*****************************************//
                for (int i = 0; i < Number_of_generators; i++)
                {
                        Elec_power[i]=Power_elc(X_vector[i].Ed_das_das,
                                                id[i],
                                                X_vector[i].Eq_das_das,
                                                iq[i],
                                                Xd2_local[i],
                                                Xq2_local[i]);
                }


                //****************************UPDATING_POWER*****************************************//
                //****************************UPDATING_POWER*****************************************//
                
                for (int i = 0; i < Number_of_generators;i++,
                                                         pointer_F_VECTOR_old++)
                {
                        F_VECTOR_CALC_single(i,
                                        pointer_F_VECTOR_old,
                                        All_data,
                                        X_vector[i],
                                        D_by_m[i],
                                        id[i],
                                        iq[i],
                                        mech_power[i],
                                        Elec_power[i]);
                }
                pointer_F_VECTOR_old = & F_vector_old[0];
                //**************************************************************************************//
                //**************************************************************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_X_VECTOR++,
                                                        pointer_F_VECTOR_old++)
                {
                        Euler_forward_wo_Exiter(pointer_to_X_elr,
                                                 pointer_X_VECTOR,
                                                 pointer_F_VECTOR_old,
                                                 del_t);
                }
                pointer_to_X_elr        = & X_elr[0];                
                pointer_X_VECTOR        = & X_vector[0];
                pointer_F_VECTOR_old    = & F_vector_old[0];

                for (int i = 0; i < Number_of_generators; 
                                                        i++,
                                                        pointer_F_VECTOR_new++)
                {
                        F_VECTOR_CALC_single(i,
                                        pointer_F_VECTOR_new,
                                        All_data,
                                        X_elr[i],
                                        D_by_m[i],
                                        id_elr[i],    // ✅ NOW USES UPDATED CURRENTS
                                        iq_elr[i],    // ✅ NOW USES UPDATED CURRENTS
                                        mech_power[i],
                                        Elec_power[i]);
                }
                pointer_F_VECTOR_new = & F_vector_new[0];
                        //f_new_calculated 
                        //888888888888888888888888888888888888888
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_to_Fsum++,
                                                        pointer_F_VECTOR_old++,
                                                        pointer_F_VECTOR_new++,
                                                        pointer_X_VECTOR++)
                {
                        Summer_diffential(      pointer_to_Fsum,
                                                pointer_F_VECTOR_old,
                                                pointer_F_VECTOR_new);
                        Sum_Xold_F_rets_Xnew
                                        (       pointer_X_VECTOR,
                                                pointer_to_Fsum,
                                                del_t);
                }
                pointer_to_Fsum         =       &F_Sum[0];
                pointer_F_VECTOR_old    =       &F_vector_old[0];
                pointer_F_VECTOR_new    =       &F_vector_new[0];
                pointer_X_VECTOR        =       &X_vector[0];

                //************** DONE TRAP 1st time *****************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_F_VECTOR_old++,
                                                        pointer_F_VECTOR_new++)
                {
                        copy_mak_set_ptr_zer(pointer_F_VECTOR_new,pointer_F_VECTOR_old);
                }
                pointer_F_VECTOR_old=&F_vector_old[0];
                pointer_F_VECTOR_new=&F_vector_new[0]; 
      
                /* duplicate reset removed - causes console spam */
        }
        //*************************************CLEARING FAULT ***************************//
        if ( t > (2*DISTURBANCES_CLEARING_TIME + FAULTED_TIME))
        {      
                // ✅ Post-fault period: Use original healthy network (no restoration needed)
                // ✅ Post-fault: Always use healthy impedance matrix
                Network_solver( All_data,
                                Number_of_generators,
                                Number_of_buses,
                                pointer_Nw_vector,
                                pointer_X_VECTOR,
                                Z_AUG_healthy,  // ✅ Use original healthy matrix
                                ptr_iq,
                                ptr_id,
                                ptr_Vt);
                //******************************
                pointer_Nw_vector=&Nw_vector[0];
                pointer_X_VECTOR=&X_vector[0];
                ptr_iq = & iq[0];
                ptr_id = & id[0];
 

                //**************************************NETWRK SOLVED ***********************************************************************************// 
                                for (int i = 0; i < Number_of_generators; i++)
                {
                        Vt[i]= sqrt (pow(Nw_vector[i].VQ,2)+pow(Nw_vector[i].VD,2));
                        vq[i]=Nw_vector[i].VQ*cos(X_vector[i].delta)+Nw_vector[i].VD*sin(X_vector[i].delta);
                        vd[i]=Nw_vector[i].VD*cos(X_vector[i].delta)-Nw_vector[i].VQ*sin(X_vector[i].delta);
                } 
                //***************************Updated Vt from new netwrok states*****************************************//
                //*****************************SOLVE FOR EFD ***************************************************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_X_VECTOR++)
                {       
                        proportional_error[i]=Vref[i]-Vt[i];
                        Efd_E_dummy_solver(   i,
                                        All_data,
                                        pointer_X_VECTOR,
                                        proportional_error[i],
                                        del_t,
                                        iq[i],
                                        Xq2_local[i],
                                        Xd2_local[i]);
                }
                pointer_X_VECTOR = & X_vector[0];
                // **************************CALCULATED EFD HERE ****************************************//
                //****************************CURRENT_UPDATED******************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        ptr_iq++,
                                                        ptr_id++)
                {
                        
                        current_update_new(X_vector[i],
                                        vd[i],
                                        vq[i],
                                        ptr_id,
                                        ptr_iq,
                                        Xd2_local[i]);
                        /* Keep E_dummy algebraically consistent with the
                           freshly updated iq */
                        X_vector[i].E_dummy = -(Xq2_local[i] - Xd2_local[i]) * iq[i];
                }
                ptr_iq=&iq[0];
                ptr_id=&id[0];
                
                //****************************CURRENT_UPDATED******************************************//
                //****************************CURRENT_UPDATED******************************************//
                
                //****************************UPDATING_POWER*****************************************//
                for (int i = 0; i < Number_of_generators; i++)
                {
                        Elec_power[i]=Power_elc(X_vector[i].Ed_das_das,
                                                id[i],
                                                X_vector[i].Eq_das_das,
                                                iq[i],
                                                Xd2_local[i],
                                                Xq2_local[i]);
                }
                //****************************UPDATING_POWER*****************************************//
                //****************************UPDATING_POWER*****************************************//
                
                for (int i = 0; i < Number_of_generators;i++,
                                                         pointer_F_VECTOR_old++)
                {
                        F_VECTOR_CALC_single(i,
                                        pointer_F_VECTOR_old,
                                        All_data,
                                        X_vector[i],
                                        D_by_m[i],
                                        id[i],
                                        iq[i],
                                        mech_power[i],
                                        Elec_power[i]);
                }
                pointer_F_VECTOR_old = & F_vector_old[0];
                //**************************************************************************************//
                //**************************************************************************************//
                for (int i = 0; i < Number_of_generators; ++i,
                                                        pointer_X_VECTOR++,
                                                        pointer_F_VECTOR_old++)
                {
                        /* 1. start from current state */
                        X_elr[i] = *pointer_X_VECTOR;

                        /* 2. predictor */
                        Euler_forward_wo_Exiter(&X_elr[i],
                                                 pointer_X_VECTOR,
                                                 pointer_F_VECTOR_old,
                                                 del_t);

                        /* 3. atomic copy-back for this generator */
                        /*
                         * DO NOT copy the Euler-predicted state back to X_vector here.
                         * Keeping X_vector unchanged until after the trapezoidal
                         * corrector preserves the two-step Heun integration scheme
                         * and prevents numerical drift/instability.
                         */
                        /* *pointer_X_VECTOR = X_elr[i]; */
                }
                
                pointer_to_X_elr=&X_elr[0];                
                pointer_X_VECTOR = & X_vector[0];
                pointer_F_VECTOR_old= &F_vector_old[0];

                for (int i = 0; i < Number_of_generators; 
                                                        i++,
                                                        pointer_F_VECTOR_new++)
                {
                        F_VECTOR_CALC_single(i,
                                        pointer_F_VECTOR_new,
                                        All_data,
                                        X_elr[i],
                                        D_by_m[i],
                                        id[i],
                                        iq[i],
                                        mech_power[i],
                                        Elec_power[i]);
                }
                pointer_F_VECTOR_new = & F_vector_new[0];
                //f_new_calculated 
                //888888888888888888888888888888888888888
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_to_Fsum++,
                                                        pointer_F_VECTOR_old++,
                                                        pointer_F_VECTOR_new++,
                                                        pointer_X_VECTOR++)
                {
                        Summer_diffential(      pointer_to_Fsum,
                                                pointer_F_VECTOR_old,
                                                pointer_F_VECTOR_new);
                        Sum_Xold_F_rets_Xnew
                                        (       pointer_X_VECTOR,
                                                pointer_to_Fsum,
                                                del_t);
                }
                pointer_to_Fsum         =       & F_Sum[0];
                pointer_F_VECTOR_old    =       &F_vector_old[0];
                pointer_F_VECTOR_new    =       &F_vector_new[0];
                pointer_X_VECTOR        =       &X_vector[0];

                fprintf(fptr,"DONE TRAP 1st time\n  ");
                fprintf(fptr," %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n ",
                                        t,
                                        X_vector[0].delta,
                                        X_vector[0].slip,
                                        X_vector[0].Efd,
                                        X_vector[0].Ed_das_das,
                                        X_vector[0].Eq_das_das,
                                        X_vector[0].Ed_das,
                                        X_vector[0].Eq_das,
                                        Nw_vector[0].VQ,
                                        Nw_vector[0].VD,
                                        vq[0],
                                        vd[0],
                                        Vt[0],
                                        X_vector[0].delta*(180/PI),
                                        proportional_error[0],
                                        iq[0],
                                        id[0],
                                        Elec_power[0]
                                        );  
                //************** DONE TRAP 1st time *****************************************//
                for (int i = 0; i < Number_of_generators; i++,
                                                        pointer_F_VECTOR_old++,
                                                        pointer_F_VECTOR_new++)
                {
                        copy_mak_set_ptr_zer(pointer_F_VECTOR_new,pointer_F_VECTOR_old);
                }
                pointer_F_VECTOR_old=&F_vector_old[0];
                pointer_F_VECTOR_new=&F_vector_new[0];
        }
        //*************************************END NORMAL OPERATION *************************//
        
        fprintf(FINAL_FILE,
                "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                                        t,
                                        X_vector[0].delta,
                                        X_vector[0].slip,
                                        X_vector[0].Efd,
                                        X_vector[0].Ed_das_das,
                                        X_vector[0].Eq_das_das,
                                        X_vector[0].Ed_das,
                                        X_vector[0].Eq_das,
                                        Nw_vector[0].VQ,
                                        Nw_vector[0].VD,
                                        Vref[0],
                                        Vt[0],
                                        proportional_error[0],
                                        iq[0],
                                        id[0],
                                        mech_power[0],
                                        Elec_power[0]
                                        );
            
    }
    printf("\n FILE_WRITTEN AS 'PLOT.csv pls SAVE IT' \n");   
}
//****************************************//
void F_VECTOR_CALC_single(int iter_num,

                        DIFF_STATES_VECTOR * pointer_to_diff_vector_iter,
                        
                        COMB_STRUCT All_data,

                        STATES X_vector_iter,

                        double D_by_m_iter,
                        double id_iter,
                        double iq_iter,
                        double P_m_iter,
                        double P_elec_iter)
{   

        double omega_base=OMEGA_BASE;

        double tdo1             =(All_data.Generator_ps[iter_num].Tdo1);
        double tqo1             =(All_data.Generator_ps[iter_num].Tqo1);
        double tdo2             =(All_data.Generator_ps[iter_num].Tdo2);
        double tqo2             =(All_data.Generator_ps[iter_num].Tqo2);
        
        double xd0_of_gen       =All_data.Generator_ps[iter_num].Xd0;
        double xq0_of_gen       =All_data.Generator_ps[iter_num].Xq0;
        double xd1_of_gen       =All_data.Generator_ps[iter_num].Xd1;
        double xq1_of_gen       =All_data.Generator_ps[iter_num].Xq1;
        double xd2_of_gen       =All_data.Generator_ps[iter_num].Xd2;
        double xq2_of_gen       =All_data.Generator_ps[iter_num].Xq2;
        double H                =All_data.Generator_ps[iter_num].H;
        double one_by_two_H     =pow(2*H,-1);// (1/(2*H));
        //MARK

        double iq           = iq_iter    ;
        double id           = id_iter    ;          

        double P_m          =   P_m_iter;
        double P_elec       =   P_elec_iter;

        double slip     =   X_vector_iter.slip;
        double Ed1      =   X_vector_iter.Ed_das;
        double Eq1      =   X_vector_iter.Eq_das;
        double Ed2      =   X_vector_iter.Ed_das_das;
        double Eq2      =   X_vector_iter.Eq_das_das;
        // double E_dummy  =   X_vector_iter.E_dummy;
        double Efd      =   X_vector_iter.Efd;
    
        //double omega_pu_elec = 1 + slip;
        double Tor_mech = P_m;
        double Tor_elec = P_elec;

        
        //define all pointer_to_diff_vectors these will be 6 for 1 generator 

        pointer_to_diff_vector_iter->f_of_Eq_das
            =(1/tdo1)*(Efd-Eq1+id*(xd0_of_gen-xd1_of_gen));

        pointer_to_diff_vector_iter->f_of_Ed_das
            =(-1/tqo1)*(Ed1  + iq * (xq0_of_gen - xq1_of_gen));

        pointer_to_diff_vector_iter->f_of_Eq_das_das
            =(1/tdo2)*( Eq1 - Eq2 + id * (xd1_of_gen - xd2_of_gen) ) ;
        
        pointer_to_diff_vector_iter->f_of_Ed_das_das
            =(1/tqo2)*(Ed1 - Ed2 + iq *(xq2_of_gen - xq1_of_gen));
        
        // 🔍 CRITICAL ROTOR DYNAMICS CALCULATION
        double torque_imbalance = Tor_mech - Tor_elec;
        double damping_term = D_by_m_iter * slip;
        pointer_to_diff_vector_iter->f_of_slip = (torque_imbalance * one_by_two_H) - damping_term;
        
        pointer_to_diff_vector_iter->f_of_delta = omega_base * slip;
}

void Efd_E_dummy_solver(int iter_num,
                        
                        COMB_STRUCT All_data,

                        STATES * X_vector,

                        double proportional_error_iter,
                        double del_t,
                        double iq_iter,
                        double Xq2_iter,
                        double Xd2_iter)
{   

        double Ka               =All_data.exciter_ps[iter_num].ka;
        double T_a              =All_data.exciter_ps[iter_num].ta;
        double del_t_by_2       = del_t*0.5;

        //state aegi read

        double iq = iq_iter;
        
        double E_dummy          =       X_vector->E_dummy; 
        double Efd              =       X_vector->Efd;

        double Efd_new ;
        double E_dummy_new ;

        double error            =       proportional_error_iter;

                double f_at_x_Edummy        =       (-ONE_By_T_DUMMY)*(E_dummy+ ((Xq2_iter-Xd2_iter)*iq) );

                double f_at_x_Efd           =       (1/T_a)*(Ka*error - Efd);
                //************EULER***********************//
                double Efd_euler        = Efd + (del_t)*(ELR_FACTOR)*(f_at_x_Efd);

                double E_dummy_euler    = E_dummy + (del_t)*(ELR_FACTOR)*(f_at_x_Edummy); 

                // sum hoga old_F , new_f mein 
        
                double f_x_plus_1_Efd = (1/T_a)*(Ka*error - Efd_euler);

                double f_x_plus_1_Edummy = (-ONE_By_T_DUMMY)*(E_dummy_euler+ ((Xq2_iter-Xd2_iter)*iq) ); ;

                double F_Sum_Efd = f_x_plus_1_Efd + f_at_x_Efd;

                double F_sum_Edummy = f_at_x_Edummy + f_x_plus_1_Edummy;

                //sum hoga X aur F_Sum
                Efd_new = Efd + del_t_by_2 *(F_Sum_Efd);

                E_dummy_new = E_dummy + del_t_by_2 * (F_sum_Edummy);

        
        //update 
        // 🔓 Remove tight ±3 pu clamp on Efd to allow the AVR to act freely
        //    (use very large limits so the clamp never activates in normal cases)
        double up_limit = 3.0;
        double low_limit = -3.0;
        
        if (Efd_new > up_limit)
        {
                Efd_new = up_limit;
        }
        if (Efd_new < low_limit)
        {
                Efd_new = low_limit;
        }
        else
        {
                Efd_new=Efd_new;
        }
        X_vector->Efd = Efd_new;
        X_vector->E_dummy=E_dummy_new;

}
//WORKS_ON_X_VECTOR
void Euler_forward_wo_Exiter(   STATES * ponter_to_X_Elr, 
                                STATES * ponter_to_X, 
                                DIFF_STATES_VECTOR * prt_F_latest,
                                double del_t_val)
{       
    double del_t = del_t_val*(ELR_FACTOR);
    ponter_to_X_Elr->Ed_das        =ponter_to_X->Ed_das       + (del_t *   prt_F_latest->f_of_Ed_das)      ;
    ponter_to_X_Elr->Eq_das        =ponter_to_X->Eq_das       + (del_t *   prt_F_latest->f_of_Eq_das)      ;
    ponter_to_X_Elr->Ed_das_das    =ponter_to_X->Ed_das_das   + (del_t *   prt_F_latest->f_of_Ed_das_das)  ;
    ponter_to_X_Elr->Eq_das_das    =ponter_to_X->Eq_das_das   + (del_t *   prt_F_latest->f_of_Eq_das_das)  ;
    ponter_to_X_Elr->delta         =ponter_to_X->delta        + (del_t *   prt_F_latest->f_of_delta)       ;
    ponter_to_X_Elr->slip          =ponter_to_X->slip         + (del_t *   prt_F_latest->f_of_slip)        ;
    ponter_to_X_Elr->E_dummy=ponter_to_X->E_dummy;
    ponter_to_X_Elr->Efd=ponter_to_X->Efd;
}
//WORKS_ON_NETWROK_VECTOR_ONLY
void Network_solver(    COMB_STRUCT DATA,
                        int Number_of_generators, 
                        int Number_of_buses ,
                        NW_STATES * pointer_nw_vector , 
                        STATES * ptr_X_vector,
                        double ** Z_AUG,
                        double * ptr_iq,
                        double * ptr_id,
                        double * ptr_Vt)
{       

        double Xd2_local[Number_of_generators];
        double Xq2_local[Number_of_generators];

        double i_Q_network[Number_of_generators];
        double i_D_network[Number_of_generators];
        
        double I_inject_NW      [2*Number_of_buses];
        double V_vec            [2*Number_of_buses];

        double iq_loc [Number_of_generators];
        // double id_loc [Number_of_generators];
       
        // ✅ FIX: Use local indices instead of modifying input pointers
        for (int i = 0; i < Number_of_generators; i++)
        {
                Xd2_local[i]=DATA.Generator_ps[i].Xd2;
                Xq2_local[i]=DATA.Generator_ps[i].Xq2;
                iq_loc[i]=ptr_iq[i]; // ✅ Use array indexing instead of pointer increment
        }

        // ✅ FIX: Use local index instead of modifying input pointer
        for (int i = 0; i < Number_of_generators; i++)
        {   

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
        for (int i = 2*Number_of_generators; i < 2*Number_of_buses; i++)
        {
                I_inject_NW[i]=0;
        }
        //I_vec_lenght 18
        double sum = 0;
        for (int i = 0; i < 2*Number_of_buses; i++)
        {       
   
                for (int j = 0; j < 2*Number_of_buses; j++)
                {
                        double prod = 0;
                        prod = Z_AUG[i][j]*I_inject_NW[j];
                        sum = sum + prod ; 
                }
                V_vec[i]=sum;
                sum = 0;
        }
        
        // ✅ FIX: Use local index instead of modifying input pointer
        for (int i = 0; i < Number_of_buses; i++)
        {
              pointer_nw_vector[i].VQ = V_vec[2*i];
              pointer_nw_vector[i].VD = V_vec[2*i + 1];  
        }
        // for (int i = 0; i < Number_of_generators; i++,
        //                                         ptr_iq++,ptr_id++)
        // {
        //         *ptr_iq=iq_loc[i];
        //         *ptr_id=id_loc[i];
        // }
        
                 
}

STATES * Sum_Xold_F_rets_Xnew(  STATES * ponter_to_X,
                                DIFF_STATES_VECTOR*pointer_to_FSUM,
                                double del_t)
{
    double del_t_by_2=(0.5)*(del_t);

    ponter_to_X->Ed_das        =ponter_to_X->Ed_das       + (del_t_by_2 *   pointer_to_FSUM->f_of_Ed_das)      ;
    ponter_to_X->Eq_das        =ponter_to_X->Eq_das       + (del_t_by_2 *   pointer_to_FSUM->f_of_Eq_das)      ;
    ponter_to_X->Ed_das_das    =ponter_to_X->Ed_das_das   + (del_t_by_2 *   pointer_to_FSUM->f_of_Ed_das_das)  ;
    ponter_to_X->Eq_das_das    =ponter_to_X->Eq_das_das   + (del_t_by_2 *   pointer_to_FSUM->f_of_Eq_das_das)  ;
    ponter_to_X->delta         =ponter_to_X->delta        + (del_t_by_2 *   pointer_to_FSUM->f_of_delta)       ;
    ponter_to_X->slip          =ponter_to_X->slip         + (del_t_by_2 *   pointer_to_FSUM->f_of_slip)        ;
    return ponter_to_X;
}


void Summer_diffential( DIFF_STATES_VECTOR * pointer_to_FSUM ,
                        DIFF_STATES_VECTOR * F_old,
                        DIFF_STATES_VECTOR * F_new)
{   
    
    pointer_to_FSUM->f_of_delta        =F_old->f_of_delta      + F_new->f_of_delta;
    pointer_to_FSUM->f_of_slip         =F_old->f_of_slip       + F_new->f_of_slip;
    pointer_to_FSUM->f_of_Efd          =F_old->f_of_Efd        + F_new->f_of_Efd;
    pointer_to_FSUM->f_of_Ed_das       =F_old->f_of_Ed_das     + F_new->f_of_Ed_das;
    pointer_to_FSUM->f_of_Eq_das       =F_old->f_of_Eq_das     + F_new->f_of_Eq_das;
    pointer_to_FSUM->f_of_Ed_das_das   =F_old->f_of_Ed_das_das + F_new->f_of_Ed_das_das;
    pointer_to_FSUM->f_of_Eq_das_das   =F_old->f_of_Eq_das_das + F_new->f_of_Eq_das_das;
    pointer_to_FSUM->f_of_E_dummy      =F_old->f_of_E_dummy    + F_new->f_of_E_dummy;

}
//set pointer used to zero
void copy_mak_set_ptr_zer(  DIFF_STATES_VECTOR * ptr_F_new,
                            DIFF_STATES_VECTOR * ptr_F_old)
{   
    
    ptr_F_old->f_of_delta        =ptr_F_new->f_of_delta      ;
    ptr_F_old->f_of_slip         =ptr_F_new->f_of_slip       ;
    ptr_F_old->f_of_Efd          =ptr_F_new->f_of_Efd        ;
    ptr_F_old->f_of_Ed_das       =ptr_F_new->f_of_Ed_das     ;
    ptr_F_old->f_of_Eq_das       =ptr_F_new->f_of_Eq_das     ;
    ptr_F_old->f_of_Ed_das_das   =ptr_F_new->f_of_Ed_das_das ;
    ptr_F_old->f_of_Eq_das_das   =ptr_F_new->f_of_Eq_das_das ;
    ptr_F_old->f_of_E_dummy      =ptr_F_new->f_of_E_dummy    ;

}


//power solver
double Torque_calc(double power , double omega)
{
    double Torque = (power / omega);
    return Torque;
}

void current_updater 
                (       
                double * ptr_iq,
                double * ptr_id,
                double vq_iter,
                double vd_iter,
                double Eq2_iter,
                double Ed2_iter,
                double Xq2_iter,
                double Xd2_iter
                )
{       

        *ptr_iq = (Ed2_iter- vd_iter )/(Xq2_iter);

        *ptr_id = (vq_iter - Eq2_iter)/(Xd2_iter);

}

void current_update_new(STATES X_iter,
                        double vd_iter,
                        double vq_iter,
                        double * ptr_id,
                        double * ptr_iq,
                        double Xd2_iter)
{
        double Eq2 = X_iter.Eq_das_das;
        double Ed2 = X_iter.Ed_das_das;
        double E_dummy = X_iter.E_dummy;
        double Ed_add_Edum = Ed2 + E_dummy;
        double vq=vq_iter;
        double vd=vd_iter;
        
        gsl_complex v_park = gsl_complex_rect (vq,vd);
        gsl_complex Z = gsl_complex_rect (0,Xd2_iter);

        gsl_complex Eq_Ed_Edum = gsl_complex_rect (Eq2,Ed_add_Edum);

        gsl_complex numerator = gsl_complex_sub (Eq_Ed_Edum,v_park);

        gsl_complex I = gsl_complex_div (numerator,Z);

        *ptr_iq=I.dat[0];
        *ptr_id=I.dat[1];

}

double modulus(double a , double b)
{       
        double mod;
        if (a>b)
        {
                mod = a-b;
        }
        if (b>a)
        {
                mod = b-a;
        }
        if (b==a)
        {
                mod = 0;
        }
        return mod;
}

void State_solver_one_machine(  int iter_num,
                                COMB_STRUCT All_data,
                                STATES * X_ptr,
                                double del_t,
                                double id_iter,
                                double iq_iter,
                                double Pm_iter,
                                double Pe_iter,
                                double D_by_M_iter,
                                double Prop_error_iter)
{              
        double omega_base=OMEGA_BASE;
        double del_t_half=(0.5)*(del_t);
        
        double tdo1             =(All_data.Generator_ps[iter_num].Tdo1);
        double tqo1             =(All_data.Generator_ps[iter_num].Tqo1);
        double tdo2             =(All_data.Generator_ps[iter_num].Tdo2);
        double tqo2             =(All_data.Generator_ps[iter_num].Tqo2);
        
        double xd0_of_gen       =All_data.Generator_ps[iter_num].Xd0;
        double xq0_of_gen       =All_data.Generator_ps[iter_num].Xq0;
        double xd1_of_gen       =All_data.Generator_ps[iter_num].Xd1;
        double xq1_of_gen       =All_data.Generator_ps[iter_num].Xq1;
        double xd2_of_gen       =All_data.Generator_ps[iter_num].Xd2;
        double xq2_of_gen       =All_data.Generator_ps[iter_num].Xq2;
        double H                =All_data.Generator_ps[iter_num].H;
        double one_by_two_H     =pow((2*H),-1);

        double Ka               =All_data.exciter_ps[iter_num].ka;
        double T_a              =All_data.exciter_ps[iter_num].ta;
        //non changin variables
        double iq           = iq_iter    ;
        double id           = id_iter    ;          

        double P_m          =   Pm_iter ;
        double P_elec       =   Pe_iter ;
        //states aegi x_0
        double delta_0=X_ptr->delta     ;
        double slip_0=X_ptr->slip       ;
        double Ed2_0=X_ptr->Ed_das_das  ;
        double Eq2_0=X_ptr->Eq_das_das  ;
        double Eq1_0=X_ptr->Eq_das      ;
        double Ed1_0=X_ptr->Ed_das      ;
        double Efd_0=X_ptr->Efd         ;

        double delta_new;
        double Ed1_new;
        double Ed2_new;
        double Efd_new;
        double Eq1_new;
        double Eq2_new;
        double slip_new;
        while   ( 
                modulus(delta_new,delta_0)>ERROR && 
                modulus(Ed1_new,Ed1_0)>ERROR && 
                modulus(Ed2_new,Ed2_0)>ERROR && 
                modulus(Efd_new,Efd_0)>ERROR && 
                modulus(Eq1_new,Eq1_0)>ERROR && 
                modulus(Eq2_new,Eq2_0)>ERROR 
                )
        {
                                
                //clc f_0 
                double omega_pu_mech_0 = 1 + slip_0;
                double Tor_mech_0 = Torque_calc(P_m,omega_pu_mech_0);
                double Tor_elec_0 = P_elec;
       
                double f_0_Eq1          =(1/tdo1)*(Efd_0-Eq1_0+id*(xd0_of_gen-xd1_of_gen))                      ;

                double f_0_Ed1          =(-1/tqo1)*(Ed1_0  + iq * (xq0_of_gen - xq1_of_gen))                    ;

                double f_0_Eq2          =(1/tdo2)*( Eq1_0 - Eq2_0 + id * (xd1_of_gen - xd2_of_gen) )            ;
        
                double f_0_Ed2          =(1/tqo2)*(Ed1_0 - Ed2_0 + iq *(xq2_of_gen - xq1_of_gen))               ;
        
                double f_0_slip         =((Tor_mech_0-Tor_elec_0)*(one_by_two_H)) - (D_by_M_iter * slip_0 )     ;
        
                double f_0_delta        =omega_base*(slip_0)                                                    ;

                double f_0_Efd          =(1/T_a)*(Ka*Prop_error_iter - Efd_0)                                   ;
                //euler hoga f_0 se X+1
                double Ed1_elr    =Ed1_0        + (     del_t * ELR_FACTOR *  f_0_Ed1   )       ;
                double Eq1_elr    =Eq1_0        + (     del_t * ELR_FACTOR *  f_0_Eq1   )       ;
                double Ed2_elr    =Ed2_0        + (     del_t * ELR_FACTOR *  f_0_Ed2   )       ;
                double Eq2_elr    =Eq2_0        + (     del_t * ELR_FACTOR *  f_0_Eq2   )       ;
                //double delta_elr  =delta_0      + (     del_t * ELR_FACTOR *  f_0_delta )       ;
                double slip_elr   =slip_0       + (     del_t * ELR_FACTOR *  f_0_slip  )       ;
                double Efd_elr    =Efd_0        + (     del_t * ELR_FACTOR *  f_0_Efd   )       ;      

                //F_new calculate hoga using ealuer wala f_+1
                double omega_pu_mech_elr = 1 + slip_elr;
                double Tor_mech_elr = Torque_calc(P_m,omega_pu_mech_elr);
                double Tor_elec_elr = P_elec;
       
                double f_1_Eq1          =(1/tdo1)*(Efd_elr-Eq1_elr+id*(xd0_of_gen-xd1_of_gen))                      ;

                double f_1_Ed1          =(-1/tqo1)*(Ed1_elr  + iq * (xq0_of_gen - xq1_of_gen))                    ;

                double f_1_Eq2          =(1/tdo2)*( Eq1_elr - Eq2_elr + id * (xd1_of_gen - xd2_of_gen) )            ;
        
                double f_1_Ed2          =(1/tqo2)*(Ed1_elr - Ed2_elr + iq *(xq2_of_gen - xq1_of_gen))               ;
        
                double f_1_slip         =((Tor_mech_elr-Tor_elec_elr)*(one_by_two_H)) - (D_by_M_iter * slip_elr )     ;
        
                double f_1_delta        =omega_base*(slip_elr)                                                    ;

                double f_1_Efd          =(1/T_a)*(Ka*Prop_error_iter - Efd_elr)                                   ;

                //f_old jo stored h wo sum hoga f_0 + f_+1
                double F_sum_delta      =f_0_delta      +       f_1_delta       ;
                double F_sum_Ed1        =f_0_Ed1        +       f_1_Ed1         ;
                double F_sum_Ed2        =f_0_Ed2        +       f_1_Ed2         ;
                double F_sum_Efd        =f_0_Efd        +       f_1_Efd         ;
                double F_sum_Eq1        =f_0_Eq1        +       f_1_Eq1         ;
                double F_sum_Eq2        =f_0_Eq2        +       f_1_Eq2         ;
                double F_sum_slip       =f_0_slip       +       f_1_slip        ;

                //X_new update hoga X2 clc karo
                delta_new       =delta_0        +del_t_half*F_sum_delta ;
                Ed1_new         =Ed1_0          +del_t_half*F_sum_Ed1   ;
                Ed2_new         =Ed2_0          +del_t_half*F_sum_Ed2   ;
                Efd_new         =Efd_0          +del_t_half*F_sum_Efd   ;
                Eq1_new         =Eq1_0          +del_t_half*F_sum_Eq1   ;
                Eq2_new         =Eq2_0          +del_t_half*F_sum_Eq2   ;
                slip_new        =slip_0         +del_t_half*F_sum_slip  ;
                
        }
        
        X_ptr->delta            = delta_new     ;
        X_ptr->Ed_das           = Ed1_new       ;
        X_ptr->Ed_das_das       =Ed2_new        ;
        X_ptr->Eq_das           =Eq1_new        ;
        X_ptr->Eq_das_das       =Eq2_new        ;
        X_ptr->slip             =slip_new       ;
}