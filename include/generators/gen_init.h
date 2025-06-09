#ifndef GEN_INIT_H
#define GEN_INIT_H

#include "../common/types_fwd.h"

/**
 * Calculate initial generator values from load flow data
 * @param NetData Combined network data structure
 * @return Pointer to allocated array of initial generator states
 */
InitialConditions * initializeGenerator(NetworkData NetData);

#endif 