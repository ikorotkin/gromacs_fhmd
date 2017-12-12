#include "data_structures.h"


double Sxyz_r(const rvec x, FHMD *fh)
{
    rvec xd;

    for(int d = 0; d < DIM; d++)
        xd[d] = x[d] - fh->box05[d];

    double r2 = xd[0]*xd[0] + xd[1]*xd[1] + xd[2]*xd[2];

    if(r2 <= fh->R12) return fh->Smin;      // S = Smin inside the radius R1
    if(r2 >= fh->R22) return fh->Smax;      // S = Smax outside the radius R2

    double r = sqrt(r2);

    return (r - fh->R1)*fh->RS + fh->Smin;                      // Linear interpolation from Smin to Smax
//  return tanh(((r - hhmd->R1)*hhmd->RS - 0.5)*6.0)*0.5 + 0.5; // Smoother interpolation
}
