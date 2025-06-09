#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "../../include/network/fault_matrix.h"
#include "../../include/common/data_types.h"

AdmittanceMatrix createFaultMatrix (int line_set, NetworkData Data) {   
    int number_of_buses = Data.constants[0].LOAD_FLOW_number;
    int number_of_lines = Data.constants[0].Number_of_lines;
    
    //***************************************************
    gsl_complex Y_bus_made [number_of_buses][number_of_buses];
    //***************************************************
    for (int i = 0; i < number_of_buses; i++) {
        for (int j = 0; j < number_of_buses; j++) {
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
    for (int i = 0; i < number_of_lines; i++) {   
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
    for (int i = 0; i < number_of_buses; i++) {
        gsl_complex SUS = gsl_complex_rect (0,(-2)*Data.trans_line_para[i].B);
        
        gsl_complex SUM =GSL_COMPLEX_ZERO;
        
        if (Data.trans_line_para[i].B == 0) {   
            //sum at ii should be zero 
            for (int j = 0; j < number_of_buses; j++) {
                SUM = gsl_complex_add (SUM,Y_bus_made[i][j]);
            }
            
            Y_bus_made[i][i]=gsl_complex_negative(SUM);
            
        }
        if (Data.trans_line_para[i].B != 0) {
            //sum at ii should b B value
            for (int j = 0; j < number_of_buses; j++) {
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
    AdmittanceMatrix Y_RTRN;
    for (int i = 0; i < number_of_buses; i++) {
        for (int j = 0; j < number_of_buses; j++) {
            Y_RTRN.MAT[i][j]=Y_bus_made[i][j];
        }
    }
    return Y_RTRN;
} 