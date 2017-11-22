#ifndef MACRO_H_
#define MACRO_H_

#include "data_structures.h"

static int I(int i, int j, int k, FH_VecI N)
{
    if(abs(i) > FHMD_IND_MAX || abs(j) > FHMD_IND_MAX || abs(k) > FHMD_IND_MAX)
    {
        printf(MAKE_RED "\nFHMD: ERROR: FH array indexes overflow: (i, j, k) = (%d, %d, %d)\n" RESET_COLOR "\n", i, j, k);
        exit(20);
    }

    while(i >= N.x) i -= N.x;
    while(j >= N.y) j -= N.y;
    while(k >= N.z) k -= N.z;
    while(i < 0)    i += N.x;
    while(j < 0)    j += N.y;
    while(k < 0)    k += N.z;

    return i + j*N.x + k*N.x*N.y;
}

#define C I(i,j,k,fh->N)            /* Current point [i][j][k] */

#endif /* MACRO_H_ */
