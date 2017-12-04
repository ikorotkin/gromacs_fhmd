#include "data_structures.h"
#include "macro.h"


void trilinear_find_neighbours(const rvec x, const int n, dvec xi, int *nbr, FHMD *fh)
{
    dvec xn;
    int  ind = fh->ind[n];
    ivec il, ir;
    dvec shift;

    PBC(xn, x, fh->box);

    /* Find 8 neighbours */
    for(int d = 0; d < DIM; d++)
    {
        if((xn[d] - fh->grid.c[ind][d]) > 0) {
            il[d]    = fh->indv[n][d];
            ir[d]    = fh->indv[n][d] + 1;
            shift[d] = -0.5*fh->grid.h[ind][d];
        } else {
            il[d]    = fh->indv[n][d] - 1;
            ir[d]    = fh->indv[n][d];
            shift[d] = 0.5*fh->grid.h[I3(fh->indv[n][0]-1, fh->indv[n][1]-1, fh->indv[n][2]-1, fh->N)][d];
        }
    }

    nbr[0] = I3(il[0], il[1], il[2], fh->N);
    nbr[1] = I3(ir[0], il[1], il[2], fh->N);
    nbr[2] = I3(il[0], ir[1], il[2], fh->N);
    nbr[3] = I3(ir[0], ir[1], il[2], fh->N);

    nbr[4] = I3(il[0], il[1], ir[2], fh->N);
    nbr[5] = I3(ir[0], il[1], ir[2], fh->N);
    nbr[6] = I3(il[0], ir[1], ir[2], fh->N);
    nbr[7] = I3(ir[0], ir[1], ir[2], fh->N);

    /* Shift and rescale coordinates */
    for(int d = 0; d < DIM; d++)
        xi[d] = (xn[d] - fh->grid.n[ind][d] + shift[d])/((fh->grid.h[nbr[0]][d] + fh->grid.h[nbr[7]][d])*0.5);
}


void trilinear_interpolation(dvec result, dvec xi, dvec f0, dvec f1, dvec f2, dvec f3, dvec f4, dvec f5, dvec f6, dvec f7)
{
    const double x = xi[0];
    const double y = xi[1];
    const double z = xi[2];

    for(int d = 0; d < DIM; d++)
    {
        result[d] = (1-z)*((1-y)*((1-x)*f0[d] + x*f1[d]) + y*((1-x)*f2[d] + x*f3[d])) +
                      (z)*((1-y)*((1-x)*f4[d] + x*f5[d]) + y*((1-x)*f6[d] + x*f7[d]));
    }
}
