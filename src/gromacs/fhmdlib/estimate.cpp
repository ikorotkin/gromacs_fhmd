#include "data_structures.h"
#include "sfunction.h"
#include "macro.h"

void fhmd_reset_statistics(FHMD *fh)
{
    MD_stat *st = &fh->stat;

    st->N            = 0;
    st->n            = 0;
    st->invN         = 0;

    st->davg_rho_md  = 0;
    st->davg_rho_fh  = 0;
    st->davg_rho2_md = 0;
    st->davg_rho2_fh = 0;

    for(int d = 0; d < DIM; d++)
    {
        st->davg_u_md[d]  = 0;
        st->davg_u2_md[d] = 0;
        st->davg_u_fh[d]  = 0;
        st->davg_u2_fh[d] = 0;
    }
}


void fhmd_collect_statistics(FHMD *fh)
{
    MD_stat   *st  = &fh->stat;
    FH_arrays *arr =  fh->arr;

    st->n++;
    double invn = 1.0/(double)(st->n);

    for(int i = 0; i < fh->Ntot; i++)
    {
        if(fh->grid.md[i] == FH_zone) continue;     // Skip pure FH region

        st->N++;

        st->avg_rho_md_cell[i] += arr[i].ro_md;
        st->avg_rho_fh_cell[i] += arr[i].ro_fh;

        st->davg_rho_md  += arr[i].ro_md;
        st->davg_rho_fh  += arr[i].ro_fh;

        // Faster convergence but less accurate
        //st->davg_rho2_md += (arr[i].ro_md - fh->total_density)*(arr[i].ro_md - fh->total_density);
        //st->davg_rho2_fh += (arr[i].ro_fh - fh->FH_dens)*(arr[i].ro_fh - fh->FH_dens);

        // More accurate estimation but much longer convergence
        st->davg_rho2_md += (arr[i].ro_md - st->avg_rho_md_cell[i]*invn)*(arr[i].ro_md - st->avg_rho_md_cell[i]*invn);
        st->davg_rho2_fh += (arr[i].ro_fh - st->avg_rho_fh_cell[i]*invn)*(arr[i].ro_fh - st->avg_rho_fh_cell[i]*invn);

        for(int d = 0; d < DIM; d++)
        {
            st->davg_u_md[d]  += arr[i].u_md[d];
            st->davg_u_fh[d]  += arr[i].u_fh[d];
            st->davg_u2_md[d] += arr[i].u_md[d]*arr[i].u_md[d];
            st->davg_u2_fh[d] += arr[i].u_fh[d]*arr[i].u_fh[d];
        }
    }

    st->invN = 1.0/(double)(st->N);
}


void fhmd_update_statistics(FHMD *fh)
{
    MD_stat *st = &fh->stat;

    st->avg_rho_md = st->davg_rho_md*st->invN;
    st->avg_rho_fh = st->davg_rho_fh*st->invN;
    st->std_rho_md = sqrt(st->davg_rho2_md*st->invN);
    st->std_rho_fh = sqrt(st->davg_rho2_fh*st->invN);

    for(int d = 0; d < DIM; d++)
    {
        st->avg_u_md[d] = st->davg_u_md[d]*st->invN;
        st->avg_u_fh[d] = st->davg_u_fh[d]*st->invN;
        st->std_u_md[d] = sqrt(fabs(st->davg_u2_md[d]*st->invN - st->avg_u_md[d]*st->avg_u_md[d]));
        st->std_u_fh[d] = sqrt(fabs(st->davg_u2_fh[d]*st->invN - st->avg_u_fh[d]*st->avg_u_fh[d]));
    }
}


void fhmd_print_statistics(FHMD *fh, t_commrec *cr)
{
    MD_stat *st = &fh->stat;

    fhmd_collect_statistics(fh);
    fhmd_update_statistics(fh);     // Every MD time step -- for stochastic integration

    if(MASTER(cr))
    {
        if((!(fh->step_MD % (fh->FH_step))    && (fh->step_MD < (fh->FH_step*10)))  ||
           (!(fh->step_MD % (fh->FH_step*10)) && (fh->step_MD < (fh->FH_step*100))) ||
            !(fh->step_MD % (fh->FH_step*100)))
        {
            if(!fh->step_MD)
            {
                printf("%8s " MAKE_LIGHT_BLUE "%10s " MAKE_BLUE "%10s " MAKE_LIGHT_BLUE "%9s %9s %9s " MAKE_BLUE "%9s %9s %9s " MAKE_LIGHT_BLUE "%9s %9s %9s "
                       MAKE_BLUE "%9s" RESET_COLOR "\n",
                       "Step", "STD_rho_MD", "STD_rho_FH", "STD_Ux_MD", "STD_Uy_MD", "STD_Uz_MD", "STD_Ux_FH", "STD_Uy_FH", "STD_Uz_FH",
                       "<Ux_MD>", "<Uy_MD>", "<Uz_MD>", "rho(C)_MD");
                printf("----------------------------------------------------------------------------------------------------------------------------------\n");
            }

            printf("\r%8d " MAKE_LIGHT_BLUE "%10.4f " MAKE_BLUE "%10.4f " MAKE_LIGHT_BLUE "%9.5f %9.5f %9.5f " MAKE_BLUE "%9.5f %9.5f %9.5f "
                   MAKE_LIGHT_BLUE "%9.2e %9.2e %9.2e" MAKE_BLUE "%10.3f",
                   fh->step_MD, st->std_rho_md, st->std_rho_fh, st->std_u_md[0], st->std_u_md[1], st->std_u_md[2],
                   st->std_u_fh[0], st->std_u_fh[1], st->std_u_fh[2], st->avg_u_md[0], st->avg_u_md[1], st->avg_u_md[2], fh->arr[fh->Ntot/2].ro_md);

            printf(RESET_COLOR "\n");

            fflush(stdout);
        }
    }
}


void fhmd_couette_avg(FHMD *fh, rvec x[], rvec v[], int N_atoms, t_commrec *cr)
{
    FH_arrays *arr = fh->arr;
    dvec       xn;
    double     S = fh->S;
    int        lr;

    const double eps = 1e-6;
    const double lrs = (double)(FHMD_COUETTE_LAYERS);

    for(int k = 0; k < FHMD_COUETTE_LAYERS; k++)
    {
        fh->avg_vel[k]   = 0;
        fh->avg_n[k]     = 0;
        fh->avg_vel_S[k] = 0;
        fh->avg_n_S[k]   = 0;
    }

    for(int n = 0; n < N_atoms; n++)
    {
        PBC(xn, x[n], fh->box);

        if(fh->S_function == moving_sphere)
            S = fhmd_Sxyz_r(x[n], fh->protein_com, fh);     // MD/FH sphere follows protein
        else if(fh->S_function == fixed_sphere)
            S = fhmd_Sxyz_r(x[n], fh->box05, fh);           // Fixed MD/FH sphere

        lr = (int)(xn[YY]/fh->box[YY]*lrs);

        if(lr < 0 || lr >= FHMD_COUETTE_LAYERS)             // This should never happen... only if the coordinates are NaN.
        {
            printf(MAKE_RED "\nFHMD: ERROR: Solution diverged. Atom #%d coordinates: (%g, %g, %g) nm\n" RESET_COLOR "\n", n, xn[0], xn[1], xn[2]);
            exit(21);
        }

        fh->avg_n[lr]++;
        fh->avg_vel[lr] += fh->vel[n][XX];

        if(S < eps)
        {
            fh->avg_n_S[lr]++;
            fh->avg_vel_S[lr] += fh->vel[n][XX];
        }
    }

    if(PAR(cr))
    {
        gmx_sumi(FHMD_COUETTE_LAYERS, fh->avg_n, cr);
        gmx_sumd(FHMD_COUETTE_LAYERS, fh->avg_vel, cr);
        gmx_sumi(FHMD_COUETTE_LAYERS, fh->avg_n_S, cr);
        gmx_sumd(FHMD_COUETTE_LAYERS, fh->avg_vel_S, cr);
    }

    for(int k = 0; k < FHMD_COUETTE_LAYERS; k++)
    {
        if(fh->avg_n[k])
            fh->avg_vel_tot[k]   += fh->avg_vel[k]/(double)(fh->avg_n[k]);
        if(fh->avg_n_S[k])
            fh->avg_vel_S_tot[k] += fh->avg_vel_S[k]/(double)(fh->avg_n_S[k]);
        fh->avg_n_tot[k]++;
    }
}

