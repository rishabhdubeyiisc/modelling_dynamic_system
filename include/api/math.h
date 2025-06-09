#ifndef MATH_API_H
#define MATH_API_H

#include "../common/common.h"

/**
 * @file math.h
 * @brief Public API for mathematical operations and linear algebra
 * 
 * This header provides interfaces for matrix operations, numerical
 * integration, and other mathematical utilities used in power system analysis.
 */

/**
 * @brief Invert a square matrix using Gauss-Jordan elimination
 * @param matrix Input matrix to invert
 * @param size Matrix dimension (size x size)
 * @param result Output: inverted matrix (must be pre-allocated)
 * @return 0 on success, non-zero on failure (singular matrix)
 * 
 * Performs in-place matrix inversion using Gauss-Jordan elimination
 * with partial pivoting for numerical stability.
 */
int invertMatrix(double** matrix, unsigned int size, double** result);

/**
 * @brief Gauss-Jordan elimination with partial pivoting
 * @param matrix Matrix to reduce (modified in-place)
 * @param rows Number of rows
 * @param cols Number of columns
 * 
 * Performs Gauss-Jordan elimination to reduce matrix to reduced
 * row echelon form. Used internally by matrix inversion routines.
 */
void gauss_jordan(double** matrix, unsigned int rows, unsigned int cols);

/**
 * @brief Swap two rows in a matrix
 * @param row1 First row index
 * @param matrix Matrix to modify
 * @param row2 Second row index
 * @param cols Number of columns
 * @return Success status
 * 
 * Utility function for row operations in matrix elimination algorithms.
 */
int swap_rows(unsigned int row1, double** matrix, unsigned int row2, unsigned int cols);

/**
 * @brief Add weighted state vectors for numerical integration
 * @param old_states Previous time step states
 * @param derivatives State derivatives
 * @param dt Time step size
 * @return Updated states pointer
 * 
 * Performs Euler forward integration: x_new = x_old + dt * dx/dt
 * Used for advancing generator state variables in time.
 */
GeneratorStates* Sum_Xold_F_rets_Xnew(
    GeneratorStates* old_states,
    StateDerivatives* derivatives,
    double dt
);

/**
 * @brief Sum differential state vectors
 * @param result Output: sum of input vectors
 * @param vector1 First input vector
 * @param vector2 Second input vector
 * 
 * Utility function for combining state derivative vectors in
 * multi-step integration schemes.
 */
void Summer_diffential(
    StateDerivatives* result,
    StateDerivatives* vector1,
    StateDerivatives* vector2
);

/**
 * @brief Copy and initialize differential state vectors
 * @param dest Destination vector
 * @param src Source vector to copy
 * 
 * Utility function for managing state derivative vectors during
 * numerical integration procedures.
 */
void copy_mak_set_ptr_zer(
    StateDerivatives* dest,
    StateDerivatives* src
);

/**
 * @brief Allocate memory for matrix inversion workspace
 * @param matrix Matrix workspace to allocate
 * @param data Network data for sizing information
 * 
 * Helper function to allocate appropriate memory for matrix
 * operations based on system size.
 */
void invert_malloc_maker(double** matrix, NetworkData data);

#endif /* MATH_API_H */ 