#include "data_structures.h"
#include "sfunction.h"
#include "macro.h"


void fhmd_update_MD_in_FH(rvec x[], rvec v[], real mass[], rvec f[], int N_atoms, FHMD *fh)
{
    FH_arrays *arr = fh->arr;
    dvec       xn;
    int        ind;
    //double     S = fh->S;

    FH_S(fh);       // Estimate S in the cells and cell faces

    /* Reset statistics */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md = 0;

        for(int d = 0; d < DIM; d++)
        {
            arr[i].uro_md[d]  = 0;
        }
    }

    /* Collect statistics */
    for(int n = 0; n < N_atoms; n++)
    {
        PBC(xn, x[n], fh->box);

        for(int d = 0; d < DIM; d++)
        {
            fh->indv[n][d] = (int)(xn[d]/fh->box[d]*(double)(fh->N[d]));
        }

        ind = I(fh->indv[n], fh->N);

        if(ind < 0 || ind >= fh->Ntot)      // This should never happen... only if the coordinates are NaN.
        {
            printf(MAKE_RED "\nFHMD: ERROR: Solution diverged. Atom #%d coordinates: (%g, %g, %g) nm\n" RESET_COLOR "\n", n, xn[0], xn[1], xn[2]);
            exit(21);
        }

        fh->ind[n] = ind;

        arr[ind].ro_md += mass[n];
/*
        if(fh->S < -1)
            S = fhmd_Sxyz_r(x[n], fh->protein_com, fh);     // MD/FH sphere follows protein
        else if(fh->S < 0)
            S = fhmd_Sxyz_r(x[n], fh->box05, fh);           // Fixed MD/FH sphere
*/
        for(int d = 0; d < DIM; d++)
        {
            arr[ind].uro_md[d] += v[n][d]*mass[n];
        }
    }

    /* Update statistics */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md *= fh->grid.ivol[i];

        for(int d = 0; d < DIM; d++)
        {
            arr[i].uro_md[d] *= fh->grid.ivol[i];
        }
    }
}


void fhmd_sum_arrays(t_commrec *cr, FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    /* Pack FHMD arrays to linear array */
    for(int i = 0; i < fh->Ntot; i++)
    {
        fh->mpi_linear[i]              = arr[i].ro_md;
        fh->mpi_linear[i + fh->Ntot]   = arr[i].uro_md[0];
        fh->mpi_linear[i + fh->Ntot*2] = arr[i].uro_md[1];
        fh->mpi_linear[i + fh->Ntot*3] = arr[i].uro_md[2];
    }

    /* Broadcast linear array */
    gmx_sumd(fh->Ntot*4, fh->mpi_linear, cr);

    /* Unpack linear array */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md     = fh->mpi_linear[i];
        arr[i].uro_md[0] = fh->mpi_linear[i + fh->Ntot];
        arr[i].uro_md[1] = fh->mpi_linear[i + fh->Ntot*2];
        arr[i].uro_md[2] = fh->mpi_linear[i + fh->Ntot*3];
    }
}


void fhmd_calculate_MDFH_terms(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    ivec ind;
    dvec alpha_term;

    for(int i = 0; i < fh->Ntot; i++)
    {
        if(arr[i].ro_md <= 0) {
            printf(MAKE_RED "\nFHMD: ERROR: Zero or NaN MD density in the cell #%d (ro_md = %g)\n" RESET_COLOR "\n", i, arr[i].ro_md);
            exit(22);
        }

        arr[i].inv_ro = 1.0/arr[i].ro_md;

        for(int d = 0; d < DIM; d++)
            arr[i].u_md[d] = arr[i].uro_md[d]*arr[i].inv_ro;

        if(fh->scheme == One_Way)
        {
            arr[i].delta_ro = arr[i].ro_fh - arr[i].ro_md;
            for(int d = 0; d < DIM; d++)
                arr[i].beta_term[d] = fh->beta*(arr[i].u_fh[d]*arr[i].ro_fh - arr[i].uro_md[d]);
        }
        else if(fh->scheme == Two_Way)
        {
//            arr[i].delta_ro = arr[i].ron_prime;
//            for(int d = 0; d < DIM; d++)
//                arr[i].beta_term[d] = fh->beta*arr[i].mn_prime[d];
            arr[i].delta_ro = arr[i].ro_fh - arr[i].ro_md;                                          // Should be layer n+1/2?
            for(int d = 0; d < DIM; d++)
                arr[i].beta_term[d] = fh->beta*(arr[i].u_fh[d]*arr[i].ro_fh - arr[i].uro_md[d]);    // Should be layer n+1/2?
        }
    }

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                {
                    arr[C].grad_ro[d] = fh->alpha*(arr[CR].delta_ro - arr[CL].delta_ro)/(0.5*(fh->grid.h[CL][d] + 2.0*fh->grid.h[C][d] + fh->grid.h[CR][d]));

                    for(int du = 0; du < DIM; du++)
                        arr[C].alpha_u_grad[du][d] = arr[C].grad_ro[d]*arr[C].S*(1 - arr[C].S)*arr[C].u_md[du];     // TODO: Fast but rough estimation!
                }
            }
        }
    }

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int du = 0; du < DIM; du++)
                {
                    for(int d = 0; d < DIM; d++)
                    {
                        alpha_term[d] = (arr[CR].alpha_u_grad[du][d] - arr[CL].alpha_u_grad[du][d])
                                /(0.5*(fh->grid.h[CL][d] + 2.0*fh->grid.h[C][d] + fh->grid.h[CR][d]));
                    }

                    arr[C].alpha_term[du] = SUM(alpha_term);
                }
            }
        }
    }
}

