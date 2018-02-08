#ifndef FHMD_FH_H_
#define FHMD_FH_H_

#include "fh_functions.h"

gmx_inline static double fmax_(double x1, double x2)
{
    return x1 > x2 ? x1 : x2;
}

gmx_inline static double fmin_(double x1, double x2)
{
    return x1 < x2 ? x1 : x2;
}

gmx_inline static double DMAX3(double x1, double x2, double x3)
{
    return fmax_(fmax_(x1,x2),x3);
}

gmx_inline static double DMIN3(double x1, double x2, double x3)
{
    return fmin_(fmin_(x1,x2),x3);
}

gmx_inline static double DRNOR()
{
    double R1  = (double)rand()/((double)(RAND_MAX)+1.0);
    double R2  = (double)rand()/((double)(RAND_MAX)+1.0);
    double R11 = sqrt(2.0*(-log(1.0 - R1)));
    double R22 = 2.0*(3.1415926535897932384626433832795)*R2;

    return R11*cos(R22);
}

void compute_random_stress(FHMD *fh);
void swap_var(FHMD *fh);


double MU, KAPPA, EOS_A, EOS_B, EOS_C, EOS_D, EOS_E, P_INIT, SOUND, VISC1, VISC2;
double T;
int    STEP;
double blend;

double T_INST, T_AVG, RHO_AVG, RHO2_AVG, P_AVG;
int    NT_AVG, N_AVG;
dvec   U_AVG, U2_AVG;

double std_rho, avg_rho, sound, avg_p, avg_T;
dvec   std_u, avg_u, avg_sound;


double HXI, HYI, HZI;

double DuxDX_CLC, DuxDX_RCC, DuxDX_RLC, DuxDX_RLB, DuxDX_RLT, DuxDX_RLD, DuxDX_RLU;
double DuyDX_CLC, DuyDX_RCC, DuyDX_RLC, DuyDX_RLB, DuyDX_RLT;
double DuzDX_CLC, DuzDX_RCC, DuzDX_RLC, DuzDX_RLD, DuzDX_RLU;

double DuxDY_TBC, DuxDY_TBL, DuxDY_TBR, DuxDY_CBC, DuxDY_TCC;
double DuyDY_TBC, DuyDY_TBL, DuyDY_TBR, DuyDY_CBC, DuyDY_TCC, DuyDY_TBD, DuyDY_TBU;
double DuzDY_CBC, DuzDY_TCC, DuzDY_TBC, DuzDY_TBD, DuzDY_TBU;

double DuxDZ_UDC, DuxDZ_UDL, DuxDZ_UDR, DuxDZ_CDC, DuxDZ_UCC;
double DuyDZ_CDC, DuyDZ_UCC, DuyDZ_UDC, DuyDZ_UDB, DuyDZ_UDT;
double DuzDZ_UDC, DuzDZ_UDL, DuzDZ_UDR, DuzDZ_UDB, DuzDZ_UDT, DuzDZ_CDC, DuzDZ_UCC;


#define HC  fh->grid.h[C][d]
#define HL  (0.5*(fh->grid.h[CL][d] + fh->grid.h[C][d]))
#define HR  (0.5*(fh->grid.h[CR][d] + fh->grid.h[C][d]))
#define SL  a[L].Sf[d]
#define SR  a[R].Sf[d]
#define SC  a[C].S


#define VISCOUS_FLUX(f) \
\
HXI = 1.0/fh->grid.h[0][0]; \
HYI = 1.0/fh->grid.h[0][1]; \
HZI = 1.0/fh->grid.h[0][2]; \
\
switch(dim) \
{ \
case 0: \
    DuxDX_RCC = (arr[I3(i+1,j  ,k  ,fh->N)].f[0] - arr[I3(i  ,j  ,k  ,fh->N)].f[0])*HXI; \
    DuxDX_CLC = (arr[I3(i  ,j  ,k  ,fh->N)].f[0] - arr[I3(i-1,j  ,k  ,fh->N)].f[0])*HXI; \
    DuyDY_TBR = (arr[I3(i+1,j+1,k  ,fh->N)].f[1] - arr[I3(i+1,j-1,k  ,fh->N)].f[1])*0.5*HYI; \
    DuyDY_TBC = (arr[I3(i  ,j+1,k  ,fh->N)].f[1] - arr[I3(i  ,j-1,k  ,fh->N)].f[1])*0.5*HYI; \
    DuyDY_TBL = (arr[I3(i-1,j+1,k  ,fh->N)].f[1] - arr[I3(i-1,j-1,k  ,fh->N)].f[1])*0.5*HYI; \
    DuzDZ_UDC = (arr[I3(i  ,j  ,k+1,fh->N)].f[2] - arr[I3(i  ,j  ,k-1,fh->N)].f[2])*0.5*HZI; \
    DuzDZ_UDR = (arr[I3(i+1,j  ,k+1,fh->N)].f[2] - arr[I3(i+1,j  ,k-1,fh->N)].f[2])*0.5*HZI; \
    DuzDZ_UDL = (arr[I3(i-1,j  ,k+1,fh->N)].f[2] - arr[I3(i-1,j  ,k-1,fh->N)].f[2])*0.5*HZI; \
\
    DuxDY_TCC = (arr[I3(i  ,j+1,k  ,fh->N)].f[0] - arr[I3(i  ,j  ,k  ,fh->N)].f[0])*HYI; \
    DuxDY_CBC = (arr[I3(i  ,j  ,k  ,fh->N)].f[0] - arr[I3(i  ,j-1,k  ,fh->N)].f[0])*HYI; \
    DuyDX_RLC = (arr[I3(i+1,j  ,k  ,fh->N)].f[1] - arr[I3(i-1,j  ,k  ,fh->N)].f[1])*0.5*HXI; \
    DuyDX_RLT = (arr[I3(i+1,j+1,k  ,fh->N)].f[1] - arr[I3(i-1,j+1,k  ,fh->N)].f[1])*0.5*HXI; \
    DuyDX_RLB = (arr[I3(i+1,j-1,k  ,fh->N)].f[1] - arr[I3(i-1,j-1,k  ,fh->N)].f[1])*0.5*HXI; \
\
    DuxDZ_CDC = (arr[I3(i  ,j  ,k  ,fh->N)].f[0] - arr[I3(i  ,j  ,k-1,fh->N)].f[0])*HZI; \
    DuxDZ_UCC = (arr[I3(i  ,j  ,k+1,fh->N)].f[0] - arr[I3(i  ,j  ,k  ,fh->N)].f[0])*HZI; \
    DuzDX_RLC = (arr[I3(i+1,j  ,k  ,fh->N)].f[2] - arr[I3(i-1,j  ,k  ,fh->N)].f[2])*0.5*HXI; \
    DuzDX_RLU = (arr[I3(i+1,j  ,k+1,fh->N)].f[2] - arr[I3(i-1,j  ,k+1,fh->N)].f[2])*0.5*HXI; \
    DuzDX_RLD = (arr[I3(i+1,j  ,k-1,fh->N)].f[2] - arr[I3(i-1,j  ,k-1,fh->N)].f[2])*0.5*HXI; \
\
    TAUR[dim][0] = MU*(VISC1*DuxDX_RCC + VISC2*(0.5*(DuyDY_TBR + DuyDY_TBC) + 0.5*(DuzDZ_UDR + DuzDZ_UDC))); \
    TAUL[dim][0] = MU*(VISC1*DuxDX_CLC + VISC2*(0.5*(DuyDY_TBC + DuyDY_TBL) + 0.5*(DuzDZ_UDC + DuzDZ_UDL))); \
    TAUR[dim][1] = MU*(      DuxDY_TCC +       (0.5*(DuyDX_RLC + DuyDX_RLT))); \
    TAUL[dim][1] = MU*(      DuxDY_CBC +       (0.5*(DuyDX_RLC + DuyDX_RLB))); \
    TAUR[dim][2] = MU*(      DuxDZ_UCC +       (0.5*(DuzDX_RLC + DuzDX_RLU))); \
    TAUL[dim][2] = MU*(      DuxDZ_CDC +       (0.5*(DuzDX_RLC + DuzDX_RLD))); \
\
    break; \
\
case 1: \
    DuyDX_RCC = (arr[I3(i+1,j  ,k  ,fh->N)].f[1] - arr[I3(i  ,j  ,k  ,fh->N)].f[1])*HXI; \
    DuyDX_CLC = (arr[I3(i  ,j  ,k  ,fh->N)].f[1] - arr[I3(i-1,j  ,k  ,fh->N)].f[1])*HXI; \
    DuxDY_TBR = (arr[I3(i+1,j+1,k  ,fh->N)].f[0] - arr[I3(i+1,j-1,k  ,fh->N)].f[0])*0.5*HYI; \
    DuxDY_TBC = (arr[I3(i  ,j+1,k  ,fh->N)].f[0] - arr[I3(i  ,j-1,k  ,fh->N)].f[0])*0.5*HYI; \
    DuxDY_TBL = (arr[I3(i-1,j+1,k  ,fh->N)].f[0] - arr[I3(i-1,j-1,k  ,fh->N)].f[0])*0.5*HYI; \
\
    DuyDY_TCC = (arr[I3(i  ,j+1,k  ,fh->N)].f[1] - arr[I3(i  ,j  ,k  ,fh->N)].f[1])*HYI; \
    DuyDY_CBC = (arr[I3(i  ,j  ,k  ,fh->N)].f[1] - arr[I3(i  ,j-1,k  ,fh->N)].f[1])*HYI ; \
    DuxDX_RLT = (arr[I3(i+1,j+1,k  ,fh->N)].f[0] - arr[I3(i-1,j+1,k  ,fh->N)].f[0])*0.5*HXI; \
    DuxDX_RLC = (arr[I3(i+1,j  ,k  ,fh->N)].f[0] - arr[I3(i-1,j  ,k  ,fh->N)].f[0])*0.5*HXI; \
    DuxDX_RLB = (arr[I3(i+1,j-1,k  ,fh->N)].f[0] - arr[I3(i-1,j-1,k  ,fh->N)].f[0])*0.5*HXI; \
    DuzDZ_UDT = (arr[I3(i  ,j+1,k+1,fh->N)].f[2] - arr[I3(i  ,j+1,k-1,fh->N)].f[2])*0.5*HZI; \
    DuzDZ_UDC = (arr[I3(i  ,j  ,k+1,fh->N)].f[2] - arr[I3(i  ,j  ,k-1,fh->N)].f[2])*0.5*HZI; \
    DuzDZ_UDB = (arr[I3(i  ,j-1,k+1,fh->N)].f[2] - arr[I3(i  ,j-1,k-1,fh->N)].f[2])*0.5*HZI; \
\
    DuyDZ_UCC = (arr[I3(i  ,j  ,k+1,fh->N)].f[1] - arr[I3(i  ,j  ,k  ,fh->N)].f[1])*HZI; \
    DuyDZ_CDC = (arr[I3(i  ,j  ,k  ,fh->N)].f[1] - arr[I3(i  ,j  ,k-1,fh->N)].f[1])*HZI; \
    DuzDY_TBU = (arr[I3(i  ,j+1,k+1,fh->N)].f[2] - arr[I3(i  ,j-1,k+1,fh->N)].f[2])*0.5*HYI; \
    DuzDY_TBC = (arr[I3(i  ,j+1,k  ,fh->N)].f[2] - arr[I3(i  ,j-1,k  ,fh->N)].f[2])*0.5*HYI; \
    DuzDY_TBD = (arr[I3(i  ,j+1,k-1,fh->N)].f[2] - arr[I3(i  ,j-1,k-1,fh->N)].f[2])*0.5*HYI; \
\
    TAUR[dim][0] = MU*(      DuyDX_RCC +       (0.5*(DuxDY_TBR + DuxDY_TBC))); \
    TAUL[dim][0] = MU*(      DuyDX_CLC +       (0.5*(DuxDY_TBL + DuxDY_TBC))); \
    TAUR[dim][1] = MU*(VISC1*DuyDY_TCC + VISC2*(0.5*(DuxDX_RLT + DuxDX_RLC) + 0.5*(DuzDZ_UDT + DuzDZ_UDC))); \
    TAUL[dim][1] = MU*(VISC1*DuyDY_CBC + VISC2*(0.5*(DuxDX_RLC + DuxDX_RLB) + 0.5*(DuzDZ_UDC + DuzDZ_UDB))); \
    TAUR[dim][2] = MU*(      DuyDZ_UCC +       (0.5*(DuzDY_TBC + DuzDY_TBU))); \
    TAUL[dim][2] = MU*(      DuyDZ_CDC +       (0.5*(DuzDY_TBC + DuzDY_TBD))); \
\
    break; \
\
case 2: \
    DuzDX_CLC = (arr[I3(i  ,j  ,k  ,fh->N)].f[2] - arr[I3(i-1,j  ,k  ,fh->N)].f[2])*HXI; \
    DuzDX_RCC = (arr[I3(i+1,j  ,k  ,fh->N)].f[2] - arr[I3(i  ,j  ,k  ,fh->N)].f[2])*HXI; \
    DuxDZ_UDL = (arr[I3(i-1,j  ,k+1,fh->N)].f[0] - arr[I3(i-1,j  ,k-1,fh->N)].f[0])*0.5*HZI; \
    DuxDZ_UDC = (arr[I3(i  ,j  ,k+1,fh->N)].f[0] - arr[I3(i  ,j  ,k-1,fh->N)].f[0])*0.5*HZI; \
    DuxDZ_UDR = (arr[I3(i+1,j  ,k+1,fh->N)].f[0] - arr[I3(i+1,j  ,k-1,fh->N)].f[0])*0.5*HZI; \
\
    DuzDY_CBC = (arr[I3(i  ,j  ,k  ,fh->N)].f[2] - arr[I3(i  ,j-1,k  ,fh->N)].f[2])*HYI; \
    DuzDY_TCC = (arr[I3(i  ,j+1,k  ,fh->N)].f[2] - arr[I3(i  ,j  ,k  ,fh->N)].f[2])*HYI; \
    DuyDZ_UDB = (arr[I3(i  ,j-1,k+1,fh->N)].f[1] - arr[I3(i  ,j-1,k-1,fh->N)].f[1])*0.5*HZI; \
    DuyDZ_UDC = (arr[I3(i  ,j  ,k+1,fh->N)].f[1] - arr[I3(i  ,j  ,k-1,fh->N)].f[1])*0.5*HZI; \
    DuyDZ_UDT = (arr[I3(i  ,j+1,k+1,fh->N)].f[1] - arr[I3(i  ,j+1,k-1,fh->N)].f[1])*0.5*HZI; \
\
    DuzDZ_CDC = (arr[I3(i  ,j  ,k  ,fh->N)].f[2] - arr[I3(i  ,j  ,k-1,fh->N)].f[2])*HZI; \
    DuzDZ_UCC = (arr[I3(i  ,j  ,k+1,fh->N)].f[2] - arr[I3(i  ,j  ,k  ,fh->N)].f[2])*HZI; \
    DuxDX_RLD = (arr[I3(i+1,j  ,k-1,fh->N)].f[0] - arr[I3(i-1,j  ,k-1,fh->N)].f[0])*0.5*HXI; \
    DuxDX_RLC = (arr[I3(i+1,j  ,k  ,fh->N)].f[0] - arr[I3(i-1,j  ,k  ,fh->N)].f[0])*0.5*HXI; \
    DuxDX_RLU = (arr[I3(i+1,j  ,k+1,fh->N)].f[0] - arr[I3(i-1,j  ,k+1,fh->N)].f[0])*0.5*HXI; \
    DuyDY_TBD = (arr[I3(i  ,j+1,k-1,fh->N)].f[1] - arr[I3(i  ,j-1,k-1,fh->N)].f[1])*0.5*HYI; \
    DuyDY_TBC = (arr[I3(i  ,j+1,k  ,fh->N)].f[1] - arr[I3(i  ,j-1,k  ,fh->N)].f[1])*0.5*HYI; \
    DuyDY_TBU = (arr[I3(i  ,j+1,k+1,fh->N)].f[1] - arr[I3(i  ,j-1,k+1,fh->N)].f[1])*0.5*HYI; \
\
    TAUR[dim][0] = MU*(      DuzDX_RCC +       (0.5*(DuxDZ_UDC + DuxDZ_UDR))); \
    TAUL[dim][0] = MU*(      DuzDX_CLC +       (0.5*(DuxDZ_UDC + DuxDZ_UDL))); \
    TAUR[dim][1] = MU*(      DuzDY_TCC +       (0.5*(DuyDZ_UDC + DuyDZ_UDT))); \
    TAUL[dim][1] = MU*(      DuzDY_CBC +       (0.5*(DuyDZ_UDC + DuyDZ_UDB))); \
    TAUR[dim][2] = MU*(VISC1*DuzDZ_UCC + VISC2*(0.5*(DuxDX_RLU + DuxDX_RLC) + 0.5*(DuyDY_TBU + DuyDY_TBC))); \
    TAUL[dim][2] = MU*(VISC1*DuzDZ_CDC + VISC2*(0.5*(DuxDX_RLC + DuxDX_RLD) + 0.5*(DuyDY_TBC + DuyDY_TBD))); \
\
    break; \
}

#endif /* FHMD_FH_H_ */
