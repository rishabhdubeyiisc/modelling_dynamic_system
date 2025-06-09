#include "../../include/common/common.h"
#include "../../include/network/y_bus_utils.h"

double ** convertComplexToReal(AdmittanceMatrix Y, NetworkData All_data) {   
    int number_of_buses = All_data.constants[0].LOAD_FLOW_number;
    
    //splitting Y_BUS
    double G_n[number_of_buses][number_of_buses];
    double B_n[number_of_buses][number_of_buses];

    for (int i = 0; i < number_of_buses; i++) {
        for (int j = 0; j < number_of_buses; j++) {
            G_n[i][j] = (Y.MAT[i][j].dat[0]);
            B_n[i][j] = (Y.MAT[i][j].dat[1]);
        }
    }
   
    Y_2X2 Y_MAT[number_of_buses][number_of_buses];
    for (int i = 0; i < number_of_buses; i++) {
        for (int j = 0; j < number_of_buses; j++) {
            Y_MAT[i][j].Y_2X2[0][0] = G_n[i][j];
            Y_MAT[i][j].Y_2X2[0][1] = (-B_n[i][j]);   
            Y_MAT[i][j].Y_2X2[1][0] = B_n[i][j];   
            Y_MAT[i][j].Y_2X2[1][1] = G_n[i][j];      
        }
    }
    
    //creating a return matrix
    double ** Y_return = (double**)malloc(2*number_of_buses*sizeof(double*));
    for (int i = 0; i < 2*number_of_buses ; i++) {    
        //allocate each row with specified columns
        Y_return[i] = (double*)malloc(2*number_of_buses*sizeof(double));
    }

    for (int i = 0; i < number_of_buses; i++) {
        for (int j = 0; j < number_of_buses; j++) {    
            Y_return[2*i][2*j]      = Y_MAT[i][j].Y_2X2[0][0];
            Y_return[2*i][2*j+1]    = Y_MAT[i][j].Y_2X2[0][1];
            Y_return[2*i+1][2*j]    = Y_MAT[i][j].Y_2X2[1][0];
            Y_return[2*i+1][2*j+1]  = Y_MAT[i][j].Y_2X2[1][1];
        }
    }
    return Y_return; 
} 