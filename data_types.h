#ifndef DATA_TYPES_H_MINE
#define DATA_TYPES_H_MINE

#include    "read_header.h"
#include    <gsl/gsl_blas.h>
#include    <gsl/gsl_complex.h>
#include    <gsl/gsl_complex_math.h>

#define PI acos(-1.0)
#define Rval 9

struct Constant_getter
{ 
  int NumberofGens,Number_of_lines,LOAD_FLOW_number;
};
typedef struct Constant_getter CONSTANT;

struct Machine_Data
{   // 0 is for Steady state,1 is for transient,2 is for Sub transient
    double  H, 
            Xd0,
            Xd1,
            Xd2,
            Xq0,
            Xq1,
            Xq2,
            Tdo1,
            Tdo2,
            Tqo1,
            Tqo2;
   
};
typedef struct Machine_Data MCD;

struct exciter
{
    double  ka, 
            ta, 
            ke, 
            te, 
            kf, 
            tf;
};
typedef struct exciter EXCT;


struct TX_PARA_all
{   
    int bus1,
        bus2;
    double 
        R,
        X,
        B;
};
typedef struct TX_PARA_all TX_PARA;


struct LOAD_FLOW_STRUCTURE
{   int bus_number;
    int bus_type;
    double  V,
            theta,
            Pg,
            Qg,
            Pl,
            Ql;
};
typedef struct LOAD_FLOW_STRUCTURE LOAD_FLOW;

struct STRUCTURE_OF_STRUCT
{   CONSTANT    * constants;
    MCD         * Generator_ps;
    EXCT        * exciter_ps;
    LOAD_FLOW   * Load_flow_ps;
    TX_PARA     * trans_line_para;
};
typedef struct STRUCTURE_OF_STRUCT COMB_STRUCT;


struct INITIAL_VALUES
{ 
    gsl_complex I_0,
                Eq_0,
                id_0,
                iq_0,
                vd_0,
                vq_0,
                Efd_0,
                Eq_dash_0,
                Ed_dash_0,
                Pe_0,
                Pm_0;
    
    double  Ed_das_das,
            Eq_das_das,
            delta_0,
            VQ,
            VD;
    
    gsl_complex I_0_DQFRAME;
};
typedef struct INITIAL_VALUES INITIAL;

struct Y_BUS_MATRIX
{
    gsl_complex MAT[Rval][Rval];
};
typedef struct Y_BUS_MATRIX Y_STRUCT;


struct STATE_VARIABLE
{
    double 
        Ed_das,
        Eq_das,
        Ed_das_das,
        Eq_das_das,
        delta,
        slip,
        Efd,
        E_dummy;
};
typedef struct STATE_VARIABLE STATES;

struct NETWORK_VARIABLES_STATES
{   
    double  VQ,
            VD;
};
typedef struct NETWORK_VARIABLES_STATES NW_STATES;

struct CURRENT_DQ
{
    double  IQ,
            ID;
};
typedef struct CURRENT_DQ I_in_DQ;


struct DIFFERENTIAL_OF_STATE 
{
    double
        f_of_Ed_das,
        f_of_Eq_das,
        f_of_Ed_das_das,
        f_of_Eq_das_das,
        f_of_delta,
        f_of_slip,
        f_of_Efd,
        f_of_E_dummy; 
};
typedef struct DIFFERENTIAL_OF_STATE DIFF_STATES_VECTOR;

struct Y_block_2_2
{
    double Y_2X2[2][2];
};
typedef struct Y_block_2_2 Y_2X2;


#endif



