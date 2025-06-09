#include "../../include/common/common.h"


NetworkData loadSystemData()
{   
    NetworkData     Array_of_pointers;
    NetworkData   * pointer_to_arr_blocks = & Array_of_pointers;
    
    SystemConstants         * constant_arr;
    GeneratorParams         * generator_arr;
    ExciterParams           * exciter_arr;
    BusData                 * Load_flow_arr;
    TransmissionLineParams  * trans_line_para_arr;
   
    
    FILE         * get_data;
    get_data     = fopen("data/system_data.txt","r");
    
    if(get_data == NULL) 
    {  
        printf("DATA NOT EXECUTED FOR MAIN DATA FILE\n");  
        exit(1); 
    };
    

    constant_arr        = malloc(sizeof(SystemConstants));
    
    fscanf(get_data,"%d",&constant_arr[0].NumberofGens);
    fscanf(get_data,"%d",&constant_arr[0].Number_of_lines);
    fscanf(get_data,"%d",&constant_arr[0].LOAD_FLOW_number);
    
    int NumberofGens        = constant_arr[0].NumberofGens;
    int Number_of_lines     = constant_arr[0].Number_of_lines;
    int LOAD_FLOW_number    = constant_arr[0].LOAD_FLOW_number;
    
    generator_arr       = malloc(NumberofGens     * sizeof(GeneratorParams));
    exciter_arr         = malloc(NumberofGens     * sizeof(ExciterParams));
    Load_flow_arr       = malloc(LOAD_FLOW_number * sizeof(BusData));
    trans_line_para_arr = malloc(Number_of_lines  * sizeof(TransmissionLineParams));
    
  
    
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

AdmittanceMatrix createAugmentedAdmittanceMatrix(NetworkData data)
{   
    FILE       * y_data;
    y_data     = fopen("data/Y_BUS.txt","r");
    
    if(y_data == NULL) 
    {  
        printf("DATA in Y BUS text NOT EXECUTED\n");  
        exit(1); 
    };
    int Bus_number=data.constants[0].LOAD_FLOW_number;
    int rows=Bus_number;
    int column=Bus_number;
    int gen_numb=data.constants[0].NumberofGens;

    AdmittanceMatrix MATRIX;
    
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

void print_matrix_complex(AdmittanceMatrix A, int order)
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
