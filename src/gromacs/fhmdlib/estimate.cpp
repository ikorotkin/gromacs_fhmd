#include "data_structures.h"


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

