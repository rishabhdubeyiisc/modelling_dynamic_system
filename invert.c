#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// HEADER FILES OF GSL
#include    <gsl/gsl_complex.h>
#include    <gsl/gsl_complex_math.h> 
// HEADER FILES I CREATED

#include "read_header.h"
#include "data_types.h"
/* */
int inv_mat(double **,unsigned int,double **);

void gauss_jordan(double **, unsigned int  ,unsigned int );

int swap_rows(unsigned int , double **, unsigned int ,unsigned int );

int inv_mat(double **inp_MAT,unsigned int order,double **op_MAT)
{
       unsigned int   i,j,k,m,n,result=0,iRow;
       double **A;
       //double **A_inverse;
       //double input_mat[2*(N+M)][2*(N+M)];
        A = (double **)malloc(order * sizeof(double *)); //Allocating rows
       	if(A == NULL) {
		printf("A memory row: allocation failed\n");
		return 99;
		}
	    for (iRow =0 ; iRow < order ; iRow++) 
        {
		    A[iRow] = (double  *)malloc(2*order * sizeof(double ));
		//Check memory validity
		    if(A[iRow] == NULL) 
            {
			printf("A memory col: allocation failed\n");
			return 99;
		}
	    }		
       n = order;
       /*
       printf("Input mat in Inverse  is\n");
       print_mat(A,n,n);
       getch();
       */
        for(i=0;i<n;i++)
       {
          for(j=0;j<2*order;j++)
          {
          	A[i][j]=0.0;
          	
          }
       }      
       for(i=0;i<n;i++)
       {
          for(j=0;j<n;j++)
          {
          	A[i][j]=inp_MAT[i][j];
          	
          }
       }
    /*--------------------------------------------------------
     * Calculating Augmented Matrix
     *-------------------------------------------------------
     */    
       for(i=0;i<n;i++)
       {
          for(j=0;j<n;j++)
          {
             if(i==j)    
                A[i][n+j]=1;
             else
                A[i][n+j]=0;
          }
       }

       gauss_jordan(A,n,2*n);


       for(i=0;i<n;i++)
       {
         // k=0;
          for(j=0;j<n;j++)
          {
             //A[i][j]=(A[i][j]/A[i][i]);
             op_MAT[i][j]=A[i][j+n];
            // k++;
          } 
       }   
       //end:
       	for (iRow =0 ; iRow < order ; iRow++) 
		 {
			free(A[iRow]);
		 }
		 free(A);
       return result;    
}


void gauss_jordan(double **A, unsigned int R ,unsigned int C)
{
	double div_val;
	unsigned int i,j,k,n=0,l=0;
	int result =0;
	
	//print_mat(A,R,R);
    for(i=0;i<R;i++)
	{
		div_val = A[i][i];
		//dbgp_printf("MAT_ARITH: Pass %d\n",i+1);
		//dbgp_printf("MAT_ARITH: div_val[%d][%d] =%f\n",i+1,i+1,div_val);
		//getch();
		if(abs(div_val)<1e-15) //if(0 == div_val)
		{
			//dbgp_printf("[MAT_ARITH :SWAPING]\n");
			result = swap_rows(i,A,R,C);
			if(1==result)
			{
				//dbgp_printf(BOLDBLUE"[MAT_ARITH: Singularity in row[%d][%d]]\n"BOLDGREEN,i+1,i+1);
				///print_mat(A,R,C);
				//system("PAUSE");
				 A[i][i]= 1e-15 ;
				 //getch();
				//exit(0);
			}
			  div_val = A[i][i];
			//dbgp_printf("[MAT_ARITH: called swap_row]\n");
		}
		for(j=0;j<C;j++)
		{
			A[i][j] =A[i][j]/div_val;
		}
			//printf("after div val G_J\n");
	        //print_mat(A,R,C);
		for(k=0;k<R;k++)
		{
			if(k!=i)
			{
				div_val =  A[k][i];
				for(j=0;j<C;j++)
				{
					A[k][j]= A[k][j] - div_val*A[i][j];
				}
			}
		}
		    //printf("after subtraction G_J\n");
	        //print_mat(A,R,C);
	}
	//gj_printf(BOLDBLUE"MAT_ARITH: GJ Completed\n");
	//#ifdef PRINT_MAT 
	//printf(BOLDGREEN"MAT_ARITH: Matrix after GJ Elimination\n");
	//print_mat(A,R,C);
	//#endif
}
/*------------------------------------------------------------------------------*/
int swap_rows(unsigned int row_num, double **A, unsigned int row,unsigned int column)
{
	unsigned int i,j;
	int result = 1;
	double* temp_row;
	temp_row = (double*) malloc( column * sizeof(temp_row));
	if (NULL != temp_row)
	{
		for(i=(row_num+1);i<row;i++)
		{
			if(0!= A[i][row_num])
			{
				result  = 0;
				for(j=row_num;j<column;j++)
				{
					temp_row[j] =A[i][j];
					A[i][j] = A[row_num][j];
					A[row_num][j]=temp_row[j];
				}
			//	printf("Row swaped b/w %d row and %d row\n",row_num+1,i+1);
				break;
			}
		}
	}
	free(temp_row);
	return(result);
}