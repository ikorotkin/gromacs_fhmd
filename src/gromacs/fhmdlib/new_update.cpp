#include "gromacs/mdtypes/group.h"
#include "gromacs/topology/atoms.h"
#include "data_structures.h"
#include "interpolation.h"


void fhmd_do_update_md(int start, int nrend,
                       double dt, int nstpcouple,
                       t_grp_tcstat *tcstat,
                       double nh_vxi[],
                       gmx_bool bNEMD, t_grp_acc *gstat, rvec accel[],
                       ivec nFreeze[],
                       real invmass[],
                       unsigned short ptype[], unsigned short cFREEZE[],
                       unsigned short cACC[], unsigned short cTC[],
                       rvec x[], rvec xprime[], rvec v[],
                       rvec f[], matrix M,
                       gmx_bool bNH, gmx_bool bPR, FHMD *fh)
{
    double imass, w_dt;
    int    gf = 0, ga = 0, gt = 0;
    rvec   vrel;
    real   vn, vv, va, vb, vnrel;
    real   lg, vxi = 0, u;
    int    n, d;

    /* FHMD variables */
    FH_arrays *arr = fh->arr;
    int        ind;
    double     invro_dt;
    double     S = fh->S;
    int        nbr[8];
    dvec       xi;
    dvec       f_fh, u_fh, alpha_term, beta_term, alpha_x_term;

    if (bNH || bPR)
    {
        /* Update with coupling to extended ensembles, used for
         * Nose-Hoover and Parrinello-Rahman coupling
         * Nose-Hoover uses the reversible leap-frog integrator from
         * Holian et al. Phys Rev E 52(3) : 2338, 1995
         */

        /* FHMD Error */
        printf(MAKE_RED "\nFHMD: ERROR: FH-MD coupling doesn't support Nose-Hoover and Parrinello-Rahman\n" RESET_COLOR "\n");
        exit(11);

    }
    else if (cFREEZE != NULL ||
             nFreeze[0][XX] || nFreeze[0][YY] || nFreeze[0][ZZ] ||
             bNEMD)
    {
        /* Update with Berendsen/v-rescale coupling and freeze or NEMD */

        /* FHMD Error */
        printf(MAKE_RED "\nFHMD: ERROR: FH-MD coupling doesn't support freeze or NEMD\n" RESET_COLOR "\n");
        exit(12);

    }
    else
    {
        /* Plain update with Berendsen/v-rescale coupling */
        for (n = start; n < nrend; n++)
        {
            if ((ptype[n] != eptVSite) && (ptype[n] != eptShell))
            {
                w_dt     = invmass[n]*dt;
                ind      = fh->ind[n];
                invro_dt = arr[ind].inv_ro*dt;

                trilinear_find_neighbours(x[n], n, xi, nbr, fh);

                trilinear_interpolation(f_fh,         xi, INTERPOLATE(f_fh));
                //trilinear_interpolation(u_fh,         xi, INTERPOLATE(u_fh));
                trilinear_interpolation(u_fh,         xi, INTERPOLATE(u_md));
                trilinear_interpolation(alpha_term,   xi, INTERPOLATE(alpha_term));
                trilinear_interpolation(beta_term,    xi, INTERPOLATE(beta_term));
                trilinear_interpolation(alpha_x_term, xi, INTERPOLATE(alpha_x_term));

#ifdef FHMD_DEBUG_INTERPOL
                if(!(n % 10000) && !(fh->step_MD % 50))
                    printf("\nStep %d, atom #%d (%g %g %g): %g %g %g\nneighbour cells: %d %d %d %d %d %d %d %d\nrescaled coordinates: %g %g %g\n",
                            fh->step_MD, n, x[n][0], x[n][1], x[n][2], u_fh[0], u_fh[1], u_fh[2],
                            nbr[0], nbr[1], nbr[2], nbr[3], nbr[4], nbr[5], nbr[6], nbr[7],
                            xi[0], xi[1], xi[2]);
#endif

                if (cTC)
                {
                    gt = cTC[n];
                }
                lg = tcstat[gt].lambda;         // Thermostat

                for (d = 0; d < DIM; d++)
                {
                    /* vn           = lg*v[n][d] + f[n][d]*w_dt; */
                    /* v[n][d]      = vn; */
                    /* xprime[n][d] = x[n][d] + vn*dt; */
                    vn           = lg*v[n][d] + ((1 - S)*f[n][d] + S*f_fh[d])*w_dt + (alpha_term[d] + beta_term[d])*invro_dt;
                    v[n][d]      = vn;
                    xprime[n][d] = x[n][d] + ((1 - S)*vn + S*u_fh[d])*dt + alpha_x_term[d]*invro_dt;
                }
            }
            else
            {
                for (d = 0; d < DIM; d++)
                {
                    v[n][d]        = 0.0;
                    xprime[n][d]   = x[n][d];
                }
            }
        }
    }
}
