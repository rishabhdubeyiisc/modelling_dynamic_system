#ifndef FAULT_MATRIX_H
#define FAULT_MATRIX_H

#include "../common/types_fwd.h"

/**
 * Create Y-bus matrix for fault conditions (line fault)
 * @param line_set Line index for fault (1-based)
 * @param Data Network data structure
 * @return Y-bus matrix with fault applied
 */
AdmittanceMatrix createFaultMatrix(int line_set, NetworkData Data);

#endif 