void partitioned_solver_copy_for_new
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
    fprintf(FINAL_FILE,"  time delta slip Efd Ed2 Eq2 Ed1 Eq1 VQ_0 VD_0 Vref Vt EROR iq_0 id_0 mech_power elec_power  \n "); 
    
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
                fprintf(fptr,"INITIAL_VALUES_TAKEN \n  ");
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
                fprintf(fptr,"NETWROK SOLVER\n  ");
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

                //**************************************NETWRK SOLVED ***********************************************************************************// 
                for (int i = 0; i < Number_of_generators; i++)
                {
                        Vt[i]= sqrt (pow(Nw_vector[i].VQ,2)+pow(Nw_vector[i].VD,2));
                        vq[i]=Nw_vector[i].VQ*cos(X_vector[i].delta)+Nw_vector[i].VD*sin(X_vector[i].delta);
                        vd[i]=Nw_vector[i].VD*cos(X_vector[i].delta)-Nw_vector[i].VQ*sin(X_vector[i].delta);
                } 
                //***************************Updated Vt from new netwrok states*****************************************//
                fprintf(fptr,"VALUES OF Vd Vq updated here\n  ");
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
                fprintf(fptr,"EFD DONE \n  ");
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
                
                fprintf(fptr,"CURRENT UPDATED \n  ");
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

                fprintf(fptr,"POWER_UPDATED \n  ");
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
                        Euler_forward_wo_Exiter(pointer_X_VECTOR,pointer_F_VECTOR_old,del_t);
                }
                                
                pointer_X_VECTOR = & X_vector[0];
                pointer_F_VECTOR_old= &F_vector_old[0];

                fprintf(fptr,"EULER DONE \n  ");
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

                for (int i = 0; i < Number_of_generators; 
                                                        i++,
                                                        pointer_F_VECTOR_new++)
                {
                        F_VECTOR_CALC_single(i,
                                        pointer_F_VECTOR_new,
                                        All_data,
                                        X_vector[i],
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
                
                initial_status = 1;
        }
        // //***********************************************************************************************************************************************************************************************************************
        // //***********************************************************************************************************************************************************************************************************************
        // // INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE
        // //***********************************************************************************************************************************************************************************************************************
        // //***********************************************************************************************************************************************************************************************************************
        // //***********************************************************************************************************************************************************************************************************************
        // // INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE INITIAL DONE
        // //***********************************************************************************************************************************************************************************************************************
        // //***********************************************************************************************************************************************************************************************************************

        //*************************************DOING DISTURBANCES *************************//
        //*************************************DOING DISTURBANCES *************************//
        //*************************************DOING DISTURBANCES *************************//
        //*************************************DOING DISTURBANCES *************************//
        //*************************************DOING DISTURBANCES *************************//
        //*************************************DOING DISTURBANCES *************************//
        //*************************************DOING DISTURBANCES *************************//        
        //*************************************DOING DISTURBANCES *************************//
        if ( t < DISTURBANCES_CLEARING_TIME || REMOVE_DISTURBANCES_FLAG == 0 )
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
                fprintf(fptr,"NETWROK SOLVER\n  ");
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

                //**************************************NETWRK SOLVED ***********************************************************************************// 
                for (int i = 0; i < Number_of_generators; i++)
                {
                        Vt[i]= sqrt (pow(Nw_vector[i].VQ,2)+pow(Nw_vector[i].VD,2));
                        vq[i]=Nw_vector[i].VQ*cos(X_vector[i].delta)+Nw_vector[i].VD*sin(X_vector[i].delta);
                        vd[i]=Nw_vector[i].VD*cos(X_vector[i].delta)-Nw_vector[i].VQ*sin(X_vector[i].delta);
                } 
                //***************************Updated Vt from new netwrok states*****************************************//
                fprintf(fptr,"VALUES OF Vd Vq updated here\n  ");
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
                fprintf(fptr,"EFD DONE \n  ");
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
                
                fprintf(fptr,"CURRENT UPDATED \n  ");
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

                fprintf(fptr,"POWER_UPDATED \n  ");
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
                        Euler_forward_wo_Exiter(pointer_X_VECTOR,pointer_F_VECTOR_old,del_t);
                }
                                
                pointer_X_VECTOR = & X_vector[0];
                pointer_F_VECTOR_old= &F_vector_old[0];

                fprintf(fptr,"EULER DONE \n  ");
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

                for (int i = 0; i < Number_of_generators; 
                                                        i++,
                                                        pointer_F_VECTOR_new++)
                {
                        F_VECTOR_CALC_single(i,
                                        pointer_F_VECTOR_new,
                                        All_data,
                                        X_vector[i],
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
                fprintf(fptr,"NETWROK SOLVER\n  ");
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

                //**************************************NETWRK SOLVED ***********************************************************************************// 
                for (int i = 0; i < Number_of_generators; i++)
                {
                        Vt[i]= sqrt (pow(Nw_vector[i].VQ,2)+pow(Nw_vector[i].VD,2));
                        vq[i]=Nw_vector[i].VQ*cos(X_vector[i].delta)+Nw_vector[i].VD*sin(X_vector[i].delta);
                        vd[i]=Nw_vector[i].VD*cos(X_vector[i].delta)-Nw_vector[i].VQ*sin(X_vector[i].delta);
                } 
                //***************************Updated Vt from new netwrok states*****************************************//
                fprintf(fptr,"VALUES OF Vd Vq updated here\n  ");
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
                fprintf(fptr,"EFD DONE \n  ");
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
                
                fprintf(fptr,"CURRENT UPDATED \n  ");
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

                fprintf(fptr,"POWER_UPDATED \n  ");
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
                        Euler_forward_wo_Exiter(pointer_X_VECTOR,pointer_F_VECTOR_old,del_t);
                }
                                
                pointer_X_VECTOR = & X_vector[0];
                pointer_F_VECTOR_old= &F_vector_old[0];

                fprintf(fptr,"EULER DONE \n  ");
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

                for (int i = 0; i < Number_of_generators; 
                                                        i++,
                                                        pointer_F_VECTOR_new++)
                {
                        F_VECTOR_CALC_single(i,
                                        pointer_F_VECTOR_new,
                                        All_data,
                                        X_vector[i],
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
        //*************************************CLEARING FAULT ***************************//
        //*************************************CLEARING FAULT ***************************//
        //*************************************CLEARING FAULT ***************************//
        //*************************************CLEARING FAULT ***************************//
        //*************************************CLEARING FAULT ***************************//
                if ( t > (2*DISTURBANCES_CLEARING_TIME + FAULTED_TIME))
        {      
                
                //***********************PASTE HERE ***********************************************//
                //***********************PASTE HERE ***********************************************//
                //***********************PASTE HERE ***********************************************//
                //***********************PASTE HERE ***********************************************//
                //***********************PASTE HERE ***********************************************//
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
                fprintf(fptr,"NETWROK SOLVER\n  ");
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

                //**************************************NETWRK SOLVED ***********************************************************************************// 
                for (int i = 0; i < Number_of_generators; i++)
                {
                        Vt[i]= sqrt (pow(Nw_vector[i].VQ,2)+pow(Nw_vector[i].VD,2));
                        vq[i]=Nw_vector[i].VQ*cos(X_vector[i].delta)+Nw_vector[i].VD*sin(X_vector[i].delta);
                        vd[i]=Nw_vector[i].VD*cos(X_vector[i].delta)-Nw_vector[i].VQ*sin(X_vector[i].delta);
                } 
                //***************************Updated Vt from new netwrok states*****************************************//
                fprintf(fptr,"VALUES OF Vd Vq updated here\n  ");
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
                fprintf(fptr,"EFD DONE \n  ");
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
                
                fprintf(fptr,"CURRENT UPDATED \n  ");
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

                fprintf(fptr,"POWER_UPDATED \n  ");
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
                        Euler_forward_wo_Exiter(pointer_X_VECTOR,pointer_F_VECTOR_old,del_t);
                }
                                
                pointer_X_VECTOR = & X_vector[0];
                pointer_F_VECTOR_old= &F_vector_old[0];

                fprintf(fptr,"EULER DONE \n  ");
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

                for (int i = 0; i < Number_of_generators; 
                                                        i++,
                                                        pointer_F_VECTOR_new++)
                {
                        F_VECTOR_CALC_single(i,
                                        pointer_F_VECTOR_new,
                                        All_data,
                                        X_vector[i],
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
        
        fprintf(FINAL_FILE,
                "% lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf   \n ",
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
    printf("\n FILE_WRITTEN AS '0_Plot.csv pls SAVE IT' \n");   
}