#include    <stdio.h>
#include    <math.h>
#include    <stdlib.h>

// HEADER FILES OF GSL
#include    <gsl/gsl_complex.h>
#include    <gsl/gsl_complex_math.h>

// HEADER FILES I CREATED

#include    "read_header.h"
#include    "data_types.h"

void next_state_single_calc(STATES * ponter_to_X, DIFF_STATES_VECTOR * prt_F_latest,double del_t)
{
    ponter_to_X->Ed_das        =ponter_to_X->Ed_das       + (del_t *   prt_F_latest->f_of_Ed_das)      ;
    ponter_to_X->Eq_das        =ponter_to_X->Eq_das       + (del_t *   prt_F_latest->f_of_Eq_das)      ;
    ponter_to_X->Ed_das_das    =ponter_to_X->Ed_das_das   + (del_t *   prt_F_latest->f_of_Ed_das_das)  ;
    ponter_to_X->Eq_das_das    =ponter_to_X->Eq_das_das   + (del_t *   prt_F_latest->f_of_Eq_das_das)  ;
    ponter_to_X->delta         =ponter_to_X->delta        + (del_t *   prt_F_latest->f_of_delta)       ;
    ponter_to_X->slip          =ponter_to_X->slip         + (del_t *   prt_F_latest->f_of_slip)        ;
    ponter_to_X->Efd           =ponter_to_X->Efd          + (del_t *   prt_F_latest->f_of_Efd)         ;
    ponter_to_X->E_dummy       =ponter_to_X->E_dummy      + (del_t *   prt_F_latest->f_of_E_dummy)     ;    
    
}

double ** Y_spitter(Y_STRUCT Y,COMB_STRUCT All_data)

{   
    int number_of_buses=All_data.constants[0].LOAD_FLOW_number;
    //splitiing Y_BUS
    double G_n[number_of_buses][number_of_buses];
    double B_n[number_of_buses][number_of_buses];

    for (int i = 0; i < number_of_buses; i++)
    {
        for (int j = 0; j < number_of_buses; j++)
        {
            G_n[i][j]=(Y.MAT[i][j].dat[0]);
            B_n[i][j]=(Y.MAT[i][j].dat[1]);
        }
        
    }
   

    Y_2X2 Y_MAT[number_of_buses][number_of_buses];
    for (int i = 0; i < number_of_buses; i++)
    {
        for (int j = 0; j < number_of_buses; j++)
        {
            Y_MAT[i][j].Y_2X2[0][0]=G_n[i][j];
            Y_MAT[i][j].Y_2X2[0][1]=(-B_n[i][j]);   
            Y_MAT[i][j].Y_2X2[1][0]=B_n[i][j];   
            Y_MAT[i][j].Y_2X2[1][1]=G_n[i][j];      
        }
        
    }
    
    
    //creating a retun matrix
   double ** Y_return=(double**)malloc(2*number_of_buses*sizeof(double*));
   for (int i = 0; i < 2*number_of_buses ; i++)
   {    //allocate each row with specified columns
       Y_return[i]=(double*)malloc(2*number_of_buses*sizeof(double));
   }

   for (int i = 0; i < number_of_buses; i++)
   {
       for (int j = 0; j < number_of_buses; j++)
        {    
            Y_return[2*i][2*j]      = Y_MAT[i][j].Y_2X2[0][0];
                
            Y_return[2*i][2*j+1]    = Y_MAT[i][j].Y_2X2[0][1];
               
            Y_return[2*i+1][2*j]    = Y_MAT[i][j].Y_2X2[1][0];
               
            Y_return[2*i+1][2*j+1]  = Y_MAT[i][j].Y_2X2[1][1];
                    
        }
       
   }
   return  Y_return; 
} 

INITIAL * Gen_data(COMB_STRUCT NetData)

{   
        
    CONSTANT   * constant    = &NetData.constants[0];
        
    int Number_of_gens       = constant[0].NumberofGens;

    INITIAL *  Gen_initial;
    Gen_initial=malloc(Number_of_gens*sizeof(INITIAL));
        
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
  
    MCD        *  Gen_array      = &NetData.Generator_ps[0];
    LOAD_FLOW  *  Load_flow_arr  = &NetData.Load_flow_ps[0];
        
    for (int x = 0; x < Number_of_gens  ; x++, Gen_array++, Load_flow_arr++)
    {
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
        // vd_0[x]     = gsl_complex_rect(( -mag_Vt_0  ) * ( sin (delta_0-theta) ),0.0) ;
        // vq_0[x]     = gsl_complex_rect((  mag_Vt_0  ) * ( cos (delta_0-theta) ),0.0) ;

        vd_0[x]     = gsl_complex_rect(( -mag_Vt_0  ) * ( sin (delta_0-angle) ),0.0) ;
        vq_0[x]     = gsl_complex_rect((  mag_Vt_0  ) * ( cos (delta_0-angle) ),0.0) ;

        Gen_initial[x].id_0=id_0[x];
        Gen_initial[x].iq_0=iq_0[x];
        Gen_initial[x].vd_0=vd_0[x];
        Gen_initial[x].vq_0=vq_0[x];


        Efd_0[x]    = gsl_complex_rect
                            (gsl_complex_abs( Eq_0[x]) - ( ( (Gen_array->Xd0  - Gen_array->Xq0 ) ) * (id_0[x].dat[0] ) )
                            ,0.0);
        //MADE CHANGES HERE THERE WAS A NEGATIVE SIGNIN BTW
        // Eq_dash_0[x]= gsl_complex_rect
        //                     (gsl_complex_abs(Efd_0[x]) - ( ( (Gen_array->Xd0) - (Gen_array->Xd1) ) * gsl_complex_abs(id_0[x] ) )
        //                     ,0.0);
        Eq_dash_0[x]= gsl_complex_rect
                            ((Efd_0[x].dat[0]) + ( ( (Gen_array->Xd0) - (Gen_array->Xd1) ) *(id_0[x].dat[0] ) )
                            ,0.0);
                    
            
        Ed_dash_0[x]= gsl_complex_rect
                            (
                                (0 - ( (  (Gen_array->Xq0) - (Gen_array->Xq1) )                    * (iq_0[x].dat[0] ) ) )
                            ,0.0);
            
        gsl_complex temp=gsl_complex_rect((Gen_array->Xd2)-(Gen_array->Xq2),0.0);

        Ed_das_das[x]=(Ed_dash_0[x].dat[0])-(iq_0[x].dat[0])*(Gen_array->Xq1-Gen_array->Xq2);
            
        Eq_das_das[x]=(Eq_dash_0[x].dat[0])+(id_0[x].dat[0])*(Gen_array->Xd1-Gen_array->Xd2);

        Pe_0[x] =gsl_complex_rect
                        ((( (((Eq_das_das[x]) ) * ( (iq_0[x].dat[0]) ) )
                        + (( (Ed_das_das[x]) ) * ((id_0[x].dat[0]) ) ) ) 
                        + (((temp.dat[0]))*((id_0[x].dat[0]))*((iq_0[x].dat[0])))),
                        0.0);
            
            
        Pm_0[x]     = Pe_0[x];

        double temp_I_real= (Eq_das_das[x]*cos(delta_0_arr[x]))-(Ed_das_das[x]*sin(delta_0_arr[x]));
        double temp_I_IMAG= (Ed_das_das[x]*cos(delta_0_arr[x]))-(Eq_das_das[x]*sin(delta_0_arr[x]));
        gsl_complex added=gsl_complex_rect(temp_I_real,temp_I_IMAG);
        gsl_complex impedance=gsl_complex_rect(0,Gen_array[x].Xd2);
        I_0_DQ_FRAME[x]= gsl_complex_div(added,impedance); 

        VD[x]= gsl_complex_abs(vd_0[x]) *cos(delta_0) - gsl_complex_abs(vq_0[x])*sin(delta_0);
        VQ[x]=  gsl_complex_abs(vq_0[x])  *cos(delta_0) + gsl_complex_abs(vd_0[x]) *sin(delta_0);  
        
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

Y_STRUCT Y_MAKER (COMB_STRUCT Data)
{   
    int number_of_buses = Data.constants[0].LOAD_FLOW_number;
    int number_of_lines = Data.constants[0].Number_of_lines;
    //***************************************************
    gsl_complex Y_bus_made [number_of_buses][number_of_buses];
    //***************************************************
    for (int i = 0; i < number_of_buses; i++)
    {
        for (int j = 0; j < number_of_buses; j++)
        {
            Y_bus_made[i][j].dat[0]=0;
            Y_bus_made[i][j].dat[1]=0;
        }    
    }
    //***************************************************
    for (int i = 0; i < number_of_lines; i++)
    {   
        int BUS1 = Data.trans_line_para[i].bus1;
        int BUS2 = Data.trans_line_para[i].bus2;

        double r = Data.trans_line_para[i].R;
        double x = Data.trans_line_para[i].X;

        gsl_complex z = gsl_complex_rect(r,x);
        gsl_complex y = gsl_complex_inverse(z); 
        gsl_complex Y = gsl_complex_negative(y);

        Y_bus_made[BUS1-1][BUS2-1]=Y;
        Y_bus_made[BUS2-1][BUS1-1]=Y;     
    }
    //*******************************************************
    for (int i = 0; i < number_of_buses; i++)
    {
        gsl_complex SUS = gsl_complex_rect (0,(-2)*Data.trans_line_para[i].B);
        
        gsl_complex SUM =GSL_COMPLEX_ZERO;
        
        if (Data.trans_line_para[i].B == 0)
        {   
            //sum at ii should be zero 
            for (int j = 0; j < number_of_buses; j++)
            {
                SUM = gsl_complex_add (SUM,Y_bus_made[i][j]);
            }
            
            Y_bus_made[i][i]=gsl_complex_negative(SUM);
            
        }
        if (Data.trans_line_para[i].B != 0)
        {
            //sum at ii should b B value
            for (int j = 0; j < number_of_buses; j++)
            {
                SUM = gsl_complex_add (SUM,Y_bus_made[i][j]);
            }
            gsl_complex neg = gsl_complex_negative(SUM);
            Y_bus_made[i][i]=gsl_complex_sub (neg,SUS);
        }
        
    }
    //***********************************************************

    // printf("\n");
    // for (int i = 0; i < 9; i++)
    // {
    //     for (int j = 0; j < 9; j++)
    //     {
    //         printf("  %lf %lfi  ",Y_bus_made[i][j].dat[0],Y_bus_made[i][j].dat[1]);
    //     }
    //     printf("\n");
    // }
    //************************************************************
    Y_STRUCT Y_RTRN;
    for (int i = 0; i < number_of_buses; i++)
    {
        for (int j = 0; j < number_of_buses; j++)
        {
            Y_RTRN.MAT[i][j]=Y_bus_made[i][j];
        }
        
    }
    return Y_RTRN;
}

Y_STRUCT Y_FAULT_MAKER (int line_set,COMB_STRUCT Data)
{   
    int number_of_buses = Data.constants[0].LOAD_FLOW_number;
    int number_of_lines = Data.constants[0].Number_of_lines;
    //***************************************************
    gsl_complex Y_bus_made [number_of_buses][number_of_buses];
    //***************************************************
    for (int i = 0; i < number_of_buses; i++)
    {
        for (int j = 0; j < number_of_buses; j++)
        {
            Y_bus_made[i][j].dat[0]=0;
            Y_bus_made[i][j].dat[1]=0;
        }    
    }
    //*******************************************************
    int BUS_FALTE_1;
    int BUS_FALTE_2;

    BUS_FALTE_1 = Data.trans_line_para[line_set-1].bus1;
    BUS_FALTE_2 = Data.trans_line_para[line_set-1].bus2;

    //***************************************************
    for (int i = 0; i < number_of_lines; i++)
    {   
        int BUS1 = Data.trans_line_para[i].bus1;
        int BUS2 = Data.trans_line_para[i].bus2;

        double r = Data.trans_line_para[i].R;
        double x = Data.trans_line_para[i].X;

        gsl_complex z = gsl_complex_rect(r,x);
        gsl_complex y = gsl_complex_inverse(z); 
        gsl_complex Y = gsl_complex_negative(y);

        Y_bus_made[BUS1-1][BUS2-1]=Y;
        Y_bus_made[BUS2-1][BUS1-1]=Y;     
    }
    //*******************************************************
    for (int i = 0; i < number_of_buses; i++)
    {
        gsl_complex SUS = gsl_complex_rect (0,(-2)*Data.trans_line_para[i].B);
        
        gsl_complex SUM =GSL_COMPLEX_ZERO;
        
        if (Data.trans_line_para[i].B == 0)
        {   
            //sum at ii should be zero 
            for (int j = 0; j < number_of_buses; j++)
            {
                SUM = gsl_complex_add (SUM,Y_bus_made[i][j]);
            }
            
            Y_bus_made[i][i]=gsl_complex_negative(SUM);
            
        }
        if (Data.trans_line_para[i].B != 0)
        {
            //sum at ii should b B value
            for (int j = 0; j < number_of_buses; j++)
            {
                SUM = gsl_complex_add (SUM,Y_bus_made[i][j]);
            }
            gsl_complex neg = gsl_complex_negative(SUM);
            Y_bus_made[i][i]=gsl_complex_sub (neg,SUS);
        }
        
    }
    //***********************************************************
    gsl_complex Y_block[2][2];
    double r_by2 = (0.5)*Data.trans_line_para[line_set-1].R;
    double x_by2 = (0.5)*Data.trans_line_para[line_set-1].X;
    double B_by2 =       Data.trans_line_para[line_set-1].B;

    gsl_complex z_by_2 = gsl_complex_rect(r_by2,x_by2);
    gsl_complex y = gsl_complex_inverse(z_by_2);
    gsl_complex B_half=gsl_complex_rect (0,B_by2);

    Y_block[0][0]=gsl_complex_add (y,B_half);
    Y_block[0][1]=Y_block[1][0]=GSL_COMPLEX_ZERO;
    Y_block[1][1]=Y_block[0][0];
    
    //*********************************************************
    Y_bus_made[BUS_FALTE_1][BUS_FALTE_1]=Y_block[0][0];
    Y_bus_made[BUS_FALTE_1][BUS_FALTE_2]=Y_block[0][1];
    Y_bus_made[BUS_FALTE_2][BUS_FALTE_1]=Y_block[1][0];
    Y_bus_made[BUS_FALTE_2][BUS_FALTE_2]=Y_block[1][1];

    //*********************************************************
    //*********************************************************
    Y_STRUCT Y_RTRN;
    for (int i = 0; i < number_of_buses; i++)
    {
        for (int j = 0; j < number_of_buses; j++)
        {
            Y_RTRN.MAT[i][j]=Y_bus_made[i][j];
        }
        
    }
    return Y_RTRN;
}

