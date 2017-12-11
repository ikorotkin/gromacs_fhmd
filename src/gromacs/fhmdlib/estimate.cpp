#include "data_structures.h"


void fhmd_reset_statistics(FHMD *fh)
{
    MD_stat *st = &fh->stat;

    for(int i = 0; i < fh->Ntot; i++)
    {
        st->N            = 0;
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
}


void fhmd_collect_statistics(FHMD *fh)
{
    MD_stat   *st  = &fh->stat;
    FH_arrays *arr =  fh->arr;

    for(int i = 0; i < fh->Ntot; i++)
    {
        st->N++;

        st->davg_rho_md  += arr[i].ro_md;
        st->davg_rho_fh  += arr[i].ro_fh;
        st->davg_rho2_md += (arr[i].ro_md - fh->total_density)*(arr[i].ro_md - fh->total_density);
        st->davg_rho2_fh += (arr[i].ro_fh - fh->total_density)*(arr[i].ro_fh - fh->total_density);

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


void fhmd_print_statistics(FHMD *fh)
{
    MD_stat *st = &fh->stat;

    fhmd_collect_statistics(fh);

    if((!(fh->step_MD %    fh->FH_step) && (fh->step_MD < 10*fh->FH_step)) ||
       (!(fh->step_MD % 10*fh->FH_step) && (fh->step_MD < 100*fh->FH_step)) ||
        !(fh->step_MD % 100*fh->FH_step))
    {
        fhmd_update_statistics(fh);

        if(!fh->step_MD)
            printf("%10s %10s %10s %9s %9s %9s %9s %9s %9s %9s %9s %9s\n",
                   "Step", "STD_rho_MD", "STD_rho_FH", "STD_Ux_MD", "STD_Uy_MD", "STD_Uz_MD", "STD_Ux_FH", "STD_Uy_FH", "STD_Uz_FH",
                   "<Ux_MD>", "<Uy_MD>", "<Uz_MD>");

        printf("\r%10d %10.4f %10.4f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.2e %9.2e %9.2e",
               fh->step_MD, st->std_rho_md, st->std_rho_fh, st->std_u_md[0], st->std_u_md[1], st->std_u_md[2],
               st->std_u_fh[0], st->std_u_fh[1], st->std_u_fh[2], st->avg_u_md[0], st->avg_u_md[1], st->avg_u_md[2]);

        printf("\n");

        fflush(stdout);
    }

}
