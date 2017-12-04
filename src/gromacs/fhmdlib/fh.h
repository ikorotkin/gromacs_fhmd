#ifndef FHMD_FH_H_
#define FHMD_FH_H_

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


void FH_init(FHMD *fh);
void FH_predictor(FHMD *fh);
void FH_corrector(FHMD *fh);
void FH_char(FHMD *fh);
void FH_do_single_timestep(FHMD *fh);
void FH_equilibrate(FHMD *fh);
void define_FH_grid(t_commrec *cr, FHMD *fh);


double MU, KAPPA, EOS_A, EOS_B, EOS_C, P_INIT, SOUND, VISC1, VISC2;
double blend;
double T_INST, T_AVG, RHO_AVG, RHO2_AVG, P_AVG;
int    NT_AVG, N_AVG;
dvec   U_AVG, U2_AVG;

#endif /* FHMD_FH_H_ */
