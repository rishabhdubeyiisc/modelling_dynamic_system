#ifndef READ_HEADER
#define READ_HEADER

#include "data_types.h"

#define OMEGA_BASE 2*PI*50
#define T_dummy 0.1
#define ONE_By_T_DUMMY 10 
#define ERROR 0.0001
#define ELR_FACTOR 0.5

;
COMB_STRUCT Reader_fn();

INITIAL * Gen_data(COMB_STRUCT);

Y_STRUCT Y_BUS(COMB_STRUCT);

double ** Y_spitter(Y_STRUCT ,COMB_STRUCT );

double ** A_MATRIX(COMB_STRUCT ,double ,int );

double ** inverse(double ,double ** );

STATES * Sum_Xold_F_rets_Xnew(STATES *,
                            DIFF_STATES_VECTOR *,
                            double );

void F_VECTOR_CALC_single(int ,

                        DIFF_STATES_VECTOR *,
                        
                        COMB_STRUCT ,

                        STATES ,

                        double ,
                        double ,
                        double ,
                        double ,
                        double );



double Power_elc(double ,double ,double ,double ,double ,double );

int inv_mat(double **,unsigned int,double **);

void gauss_jordan(double **, unsigned int  ,unsigned int );

int swap_rows(unsigned int , double **, unsigned int ,unsigned int );

void invert_malloc_maker(double ** ,COMB_STRUCT);

void Summer_diffential( DIFF_STATES_VECTOR * ,
                        DIFF_STATES_VECTOR * ,
                        DIFF_STATES_VECTOR * );

Y_STRUCT Y_BUS_AUG(COMB_STRUCT );

double ** Z_splitter(COMB_STRUCT ,char []);

void print_matrix_octave(double ** , int );

void next_state_single_calc(STATES * , DIFF_STATES_VECTOR * ,double );

double Torque_calc(double , double);

void Euler_forward_wo_Exiter(STATES*,STATES * , DIFF_STATES_VECTOR * ,double );

void Network_solver(    COMB_STRUCT ,
                        int , 
                        int  ,
                        NW_STATES * , 
                        STATES * ,
                        double **,
                        double * ,
                        double * ,
                        double *);

void partitioned_solver_new
        (
        COMB_STRUCT ,
        INITIAL * ,
        double ,
        double ,
        double ,
        double ,
        int ,
        int,
        FILE *,
        double ** ,
        double ** ,
        Y_STRUCT ,
        FILE * ,
        double,
        double,
        double,
        int,
        double
        );

void current_updater 
                (       
                double *,
                double *,
                double ,
                double ,
                double ,
                double ,
                double ,
                double 
                );
void copy_mak_set_ptr_zer(  DIFF_STATES_VECTOR * ,
                            DIFF_STATES_VECTOR * );

void print_matrix_complex(Y_STRUCT , int );

Y_STRUCT Y_MAKER (COMB_STRUCT );

Y_STRUCT Y_BUS_AUG_MAKER_NOT_FROM_DATA(Y_STRUCT ,COMB_STRUCT);
Y_STRUCT Y_FAULT_MAKER (int ,COMB_STRUCT );

void Efd_E_dummy_solver(int ,
                        
                        COMB_STRUCT ,

                        STATES * ,

                        double ,
                        double ,
                        double ,
                        double ,
                        double );

void current_update_new(STATES ,
                        double,
                        double ,
                        double * ,
                        double * ,
                        double );

double modulus(double , double );

void State_solver_one_machine(  int ,
                                COMB_STRUCT ,
                                STATES * ,
                                double ,
                                double ,
                                double,
                                double ,
                                double,
                                double ,
                                double );


#endif
