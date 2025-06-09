#ifndef Y_BUS_UTILS_H
#define Y_BUS_UTILS_H

#include "../common/types_fwd.h"

/**
 * Split complex Y-bus matrix into real 2x2 block format
 * @param Y Complex Y-bus matrix
 * @param All_data Network data structure
 * @return Pointer to allocated real matrix in 2x2 block format
 */
double ** convertComplexToReal(AdmittanceMatrix Y, NetworkData All_data);

#endif 