#include    <stdio.h>
#include    <math.h>
#include    <stdlib.h>

// HEADER FILES I CREATED

#include    "read_header.h"
#include    "data_types.h"


COMB_STRUCT Reader_fn()
{   
    COMB_STRUCT     Array_of_pointers;
    COMB_STRUCT   * pointer_to_arr_blocks = & Array_of_pointers;
    
    CONSTANT     * constant_arr;
    MCD          * generator_arr;
    EXCT         * exciter_arr;
    LOAD_FLOW    * Load_flow_arr;
    TX_PARA      * trans_line_para_arr;
   
    
    FILE         * get_data;
    get_data     = fopen("RISHABH_DATA.txt","r");
    
    if(get_data == NULL) 
    {  
        printf("DATA NOT EXECUTED FOR MAIN DATA FILE\n");  
        exit(1); 
    };
    

    constant_arr        = malloc(sizeof(CONSTANT));
    
    fscanf(get_data,"%d",&constant_arr[0].NumberofGens);
    fscanf(get_data,"%d",&constant_arr[0].Number_of_lines);
    fscanf(get_data,"%d",&constant_arr[0].LOAD_FLOW_number);
    
    int NumberofGens        = constant_arr[0].NumberofGens;
    int Number_of_lines     = constant_arr[0].Number_of_lines;
    int LOAD_FLOW_number    = constant_arr[0].LOAD_FLOW_number;
    
    generator_arr       = malloc(NumberofGens     * sizeof(MCD));
    exciter_arr         = malloc(NumberofGens     * sizeof(EXCT));
    Load_flow_arr       = malloc(LOAD_FLOW_number * sizeof(LOAD_FLOW));
    trans_line_para_arr = malloc(Number_of_lines  * sizeof(TX_PARA));
    
  
    
    for ( int i = 0; i < NumberofGens ; i++)
    {   
        fscanf(get_data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf" ,
                &generator_arr[i].H,
                &generator_arr[i].Xd0,
                &generator_arr[i].Xd1,
                &generator_arr[i].Xd2,
                &generator_arr[i].Xq0,
                &generator_arr[i].Xq1,
                &generator_arr[i].Xq2,
                &generator_arr[i].Tdo1,
                &generator_arr[i].Tdo2,
                &generator_arr[i].Tqo1,
                &generator_arr[i].Tqo2);
    
    }
    
    for (int i = 0; i < NumberofGens; i++)
    {
        fscanf( get_data,"%lf %lf %lf %lf %lf %lf",
                &exciter_arr[i].ka,
                &exciter_arr[i].ta,
                &exciter_arr[i].ke,
                &exciter_arr[i].te,
                &exciter_arr[i].kf,
                &exciter_arr[i].tf 
                                    );
    }
    
    for (int i = 0; i < Number_of_lines; i++)
    {
         fscanf(get_data,"%d %d %lf %lf %lf",
                &trans_line_para_arr[i].bus1,
                &trans_line_para_arr[i].bus2,
                &trans_line_para_arr[i].R,
                &trans_line_para_arr[i].X,
                &trans_line_para_arr[i].B
                                            );
    }
    
    for (int i = 0; i < LOAD_FLOW_number; i++)
    {
        fscanf( get_data,"%d %d %lf %lf %lf %lf %lf %lf",
                &Load_flow_arr[i].bus_number,
                &Load_flow_arr[i].bus_type,
                &Load_flow_arr[i].V,
                &Load_flow_arr[i].theta,
                &Load_flow_arr[i].Pg,
                &Load_flow_arr[i].Qg,
                &Load_flow_arr[i].Pl,
                &Load_flow_arr[i].Ql
                                    );
    }
    
   
    pointer_to_arr_blocks->constants        = constant_arr         ;
    pointer_to_arr_blocks->Generator_ps     = generator_arr        ;
    pointer_to_arr_blocks->exciter_ps       = exciter_arr          ;
    pointer_to_arr_blocks->Load_flow_ps     = Load_flow_arr        ;
    pointer_to_arr_blocks->trans_line_para  = trans_line_para_arr  ;

    //new here

   
    fclose(get_data)        ;
    printf("\n EXECUTED READING FROM FILE \n");
    return Array_of_pointers;
};

Y_STRUCT Y_BUS(COMB_STRUCT data)
{   
    FILE       * y_data;
    y_data     = fopen("Y_BUS.txt","r");
    
    if(y_data == NULL) 
    {  
        printf("DATA in Y BUS text NOT EXECUTED\n");  
        exit(1); 
    };
    int Bus_number=data.constants[0].LOAD_FLOW_number;
    int rows=Bus_number;
    int column=Bus_number;

    Y_STRUCT MATRIX;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < column; j++)
        {
            fscanf(y_data,"%lf %lf",&MATRIX.MAT[i][j].dat[0],&MATRIX.MAT[i][j].dat[1]);
        }
        
    }
    
    printf("EXECUTED MAKING Y_BUS from file");
    return MATRIX;
    

};

Y_STRUCT Y_BUS_AUG(COMB_STRUCT data)
{   
    FILE       * y_data;
    y_data     = fopen("Y_BUS.txt","r");
    
    if(y_data == NULL) 
    {  
        printf("DATA in Y BUS text NOT EXECUTED\n");  
        exit(1); 
    };
    int Bus_number=data.constants[0].LOAD_FLOW_number;
    int rows=Bus_number;
    int column=Bus_number;
    int gen_numb=data.constants[0].NumberofGens;

    Y_STRUCT MATRIX;
    
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < column; j++)
        {
            fscanf(y_data,"%lf %lf",&MATRIX.MAT[i][j].dat[0],&MATRIX.MAT[i][j].dat[1]);
        }
        
    }

    //modify Y BUS
    for (int i = 0; i < gen_numb; i++)
    {   
        double reactance = data.Generator_ps[i].Xd2;
        gsl_complex impedance   = gsl_complex_rect (0,reactance);
        gsl_complex one         = gsl_complex_rect (1,0);
        gsl_complex admittance  = gsl_complex_div(one,impedance);
        
        MATRIX.MAT[i][i]=gsl_complex_add(MATRIX.MAT[i][i],admittance);
    }
    
    for (int i = gen_numb; i < rows; i++)
    {
        gsl_complex power       =   gsl_complex_rect(data.Load_flow_ps[i].Pl,0) ;
        gsl_complex reactive    =   gsl_complex_rect(0,data.Load_flow_ps[i].Ql) ;
        double vt_sqr_real      =   pow(data.Load_flow_ps[i].V,2);
        gsl_complex vt_sqr      =   gsl_complex_rect(vt_sqr_real,0) ;

        gsl_complex P_minus_jQ  =   gsl_complex_sub (power,reactive) ;

        gsl_complex load_admit  =   gsl_complex_div(P_minus_jQ,vt_sqr);

        MATRIX.MAT[i][i]        =   gsl_complex_add (MATRIX.MAT[i][i],load_admit);
    }
    
    
    printf("EXECUTED MAKING Y_BUS from file");
    return MATRIX;
    
};

double ** Z_splitter(COMB_STRUCT data,char file_name[])
{   
    double ** Z_split=(double**)malloc(2*data.constants[0].LOAD_FLOW_number *sizeof(double*));
    for (int i = 0; i < 2*data.constants[0].LOAD_FLOW_number ; i++)
    { 
       Z_split[i]=(double*)malloc(2*data.constants[0].LOAD_FLOW_number*sizeof(double));
    }
    
    FILE *Z_sp_ptr;
    Z_sp_ptr = fopen(file_name,"r");
    if(Z_sp_ptr == NULL)
    {
      printf("Z_____inverted_mat.txt Error!");   
      exit(1);             
    }
    for (int i = 0; i < 2*data.constants[0].LOAD_FLOW_number; i++)
    {
      for (int j = 0; j < 2*data.constants[0].LOAD_FLOW_number; j++)
      { 
          fscanf(Z_sp_ptr," %lf ",&Z_split[i][j]);
      }
    }
    fclose(Z_sp_ptr);
    return Z_split;

}

void print_matrix_octave(double ** matrix, int order)
{
    FILE *file;

    file = fopen("print_mat_octave.txt","w");
    if(file == NULL)
    {
      printf("Error!");   
      exit(1);             
    } 
    
    for (int i = 0; i < order; i++)
    {
        for (int j = 0; j < order; j++)
        {   
            if ( j == order -1)
            {
               fprintf(file," %lf  ", matrix[i][j]); 
            }
            else
            {
                fprintf(file," %lf , ", matrix[i][j]);
            }
            
        }
        fprintf(file," ;\n ");
    }
    fclose(file);

}

void print_matrix_complex(Y_STRUCT A, int order)
{
    FILE *file;

    file = fopen("print_complex_mat.txt","w");
    if(file == NULL)
    {
      printf("Error!");   
      exit(1);             
    } 
    
    for (int i = 0; i < order; i++)
    {
        for (int j = 0; j < order; j++)
        {   
            if ( j == order -1)
            {
               fprintf(file," %lf %lfi  ", A.MAT[i][j].dat[0],A.MAT[i][j].dat[1]); 
            }
            else
            {
                fprintf(file," %lf %lfi , ",A.MAT[i][j].dat[0],A.MAT[i][j].dat[1]);
            }
            
        }
        fprintf(file," ;\n ");
    }
    fclose(file);

}

Y_STRUCT Y_BUS_AUG_MAKER_NOT_FROM_DATA(Y_STRUCT Y_BUS,COMB_STRUCT data)
{   

    int Bus_number=data.constants[0].LOAD_FLOW_number;
    int rows=Bus_number;
    int column=Bus_number;
    int gen_numb=data.constants[0].NumberofGens;

    Y_STRUCT MATRIX;
    
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < column; j++)
        {
            MATRIX.MAT[i][j]=Y_BUS.MAT[i][j];
        }
        
    }

    //modify Y BUS
    for (int i = 0; i < gen_numb; i++)
    {   
        double reactance = data.Generator_ps[i].Xd2;
        gsl_complex impedance   = gsl_complex_rect (0,reactance);
        gsl_complex one         = gsl_complex_rect (1,0);
        gsl_complex admittance  = gsl_complex_div(one,impedance);
        
        MATRIX.MAT[i][i]=gsl_complex_add(MATRIX.MAT[i][i],admittance);
    }
    
    for (int i = gen_numb; i < rows; i++)
    {
        gsl_complex power       =   gsl_complex_rect(data.Load_flow_ps[i].Pl,0) ;
        gsl_complex reactive    =   gsl_complex_rect(0,data.Load_flow_ps[i].Ql) ;
        double vt_sqr_real      =   pow(data.Load_flow_ps[i].V,2);
        gsl_complex vt_sqr      =   gsl_complex_rect(vt_sqr_real,0) ;

        gsl_complex P_minus_jQ  =   gsl_complex_sub (power,reactive) ;

        gsl_complex load_admit  =   gsl_complex_div(P_minus_jQ,vt_sqr);

        MATRIX.MAT[i][i]        =   gsl_complex_add (MATRIX.MAT[i][i],load_admit);
    }
    
    
    printf("\nEXECUTED MAKING Y_BUS MADE FROM TRANSMISSION LINE DATA\n");
    return MATRIX;
    
};
