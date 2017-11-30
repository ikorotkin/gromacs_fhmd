#ifndef MACRO_H_
#define MACRO_H_

#include "data_structures.h"

#define NX      fh->N[0]                /* Number of FH cells along X axis */
#define NY      fh->N[1]                /* Number of FH cells along Y axis */
#define NZ      fh->N[2]                /* Number of FH cells along Z axis */

#define C       I(ind, fh->N)           /* Point [i][j][k] */


#define ASSIGN_IND(ind, i, j, k) \
    ind[0] = i; \
    ind[1] = j; \
    ind[2] = k;


static void PBC(dvec xn, const rvec x, const dvec box)
{
    for(int d = 0; d < DIM; d++)
    {
        xn[d] = x[d];

        if(fabs(xn[d]) > FHMD_MAX_LENGTH)
        {
            printf(MAKE_RED "\nFHMD: ERROR: Solution diverged. Atom's coordinates: (%g, %g, %g) nm\n" RESET_COLOR "\n", x[0], x[1], x[2]);
            exit(20);
        }

        while(xn[d] < 0)       xn[d] += box[d];
        while(xn[d] >= box[d]) xn[d] -= box[d];
    }
}


static int I(const ivec ind, const ivec N)
{
    ivec indn;

    copy_ivec(ind, indn);

    for(int d = 0; d < DIM; d++)
    {
        if(indn[d] < 0)     indn[d] += N[d];
        if(indn[d] >= N[d]) indn[d] -= N[d];
    }

    return indn[0] + indn[1]*N[0] + indn[2]*N[0]*N[1];
}


static int I3(const int i, const int j, const int k, const ivec N)
{
    ivec ind;

    ASSIGN_IND(ind, i, j, k);

    return I(ind, N);
}

#endif /* MACRO_H_ */
