#include "data_structures.h"
#include "sfunction.h"
#include "macro.h"


void fhmd_update_MD_in_FH(rvec x[], rvec v[], real mass[], rvec f[], int N_atoms, FHMD *fh)
{
    FH_arrays *arr = fh->arr;
    dvec       xn;
    int        ind;
    double     S = fh->S;

    /* Reset statistics */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md   = 0;
        arr[i].ro_md_s = 0;

        for(int d = 0; d < DIM; d++)
        {
            arr[i].uro_md[d]    = 0;
            arr[i].uro_md_s[d]  = 0;
        }
    }

    /* Collect statistics */
    for(int n = 0; n < N_atoms; n++)
    {
        PBC(xn, x[n], fh->box);

        for(int d = 0; d < DIM; d++)
        {
            fh->indv[n][d] = (int)(xn[d]/fh->box[d]*(double)(fh->N_md[d])) + fh->N_shift[d];
        }

        ind = I(fh->indv[n], fh->N);

        if(ind < 0 || ind >= fh->Ntot)      // This should never happen... only if the coordinates are NaN.
        {
            printf(MAKE_RED "\nFHMD: ERROR: Solution diverged. Atom #%d coordinates: (%g, %g, %g) nm\n" RESET_COLOR "\n", n, xn[0], xn[1], xn[2]);
            exit(21);
        }

        fh->ind[n] = ind;

        if(fh->S_function == moving_sphere)
            S = fhmd_Sxyz_r(x[n], fh->protein_com, fh);     // MD/FH sphere follows protein
        else if(fh->S_function == fixed_sphere)
            S = fhmd_Sxyz_r(x[n], fh->box05, fh);           // Fixed MD/FH sphere

        arr[ind].ro_md   += mass[n];
        arr[ind].ro_md_s += (1 - S)*mass[n];

        for(int d = 0; d < DIM; d++)
        {
            arr[ind].uro_md[d]   += v[n][d]*mass[n];
            arr[ind].uro_md_s[d] += (1 - S)*v[n][d]*mass[n];
        }
    }

    /* Update statistics */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md   *= fh->grid.ivol[i];
        arr[i].ro_md_s *= fh->grid.ivol[i];

        for(int d = 0; d < DIM; d++)
        {
            arr[i].uro_md[d]   *= fh->grid.ivol[i];
            arr[i].uro_md_s[d] *= fh->grid.ivol[i];
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

        fh->mpi_linear[i + fh->Ntot*4] = arr[i].ro_md_s;
        fh->mpi_linear[i + fh->Ntot*5] = arr[i].uro_md_s[0];
        fh->mpi_linear[i + fh->Ntot*6] = arr[i].uro_md_s[1];
        fh->mpi_linear[i + fh->Ntot*7] = arr[i].uro_md_s[2];
    }

    /* Broadcast linear array */
    gmx_sumd(fh->Ntot*8, fh->mpi_linear, cr);

    /* Unpack linear array */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md     = fh->mpi_linear[i];
        arr[i].uro_md[0] = fh->mpi_linear[i + fh->Ntot];
        arr[i].uro_md[1] = fh->mpi_linear[i + fh->Ntot*2];
        arr[i].uro_md[2] = fh->mpi_linear[i + fh->Ntot*3];

        arr[i].ro_md_s     = fh->mpi_linear[i + fh->Ntot*4];
        arr[i].uro_md_s[0] = fh->mpi_linear[i + fh->Ntot*5];
        arr[i].uro_md_s[1] = fh->mpi_linear[i + fh->Ntot*6];
        arr[i].uro_md_s[2] = fh->mpi_linear[i + fh->Ntot*7];
    }
}


void fhmd_calculate_MDFH_terms(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    ivec ind;
    dvec alpha_term;
    int  ij;

    for(int k = fh->N_shift[2]; k < (fh->N_md[2] + fh->N_shift[2]); k++)
    {
        for(int j = fh->N_shift[1]; j < (fh->N_md[1] + fh->N_shift[1]); j++)
        {
            for(int i = fh->N_shift[0]; i < (fh->N_md[0] + fh->N_shift[0]); i++)
            {
                ASSIGN_IND(ind, i, j, k);

                if(arr[C].ro_md <= 0) {
                    printf(MAKE_RED "\nFHMD: ERROR: Zero or NaN MD density in the cell %d-%d-%d (ro_md = %g)\n" RESET_COLOR "\n", i, j, k, arr[C].ro_md);
                    exit(22);
                }

                arr[C].inv_ro = 1.0/arr[C].ro_md;

                for(int d = 0; d < DIM; d++)
                    arr[C].u_md[d] = arr[C].uro_md[d]*arr[C].inv_ro;

                if(fh->scheme == One_Way)
                {
                    arr[C].delta_ro = arr[C].ro_fh - arr[C].ro_md;
                    for(int d = 0; d < DIM; d++)
                        arr[C].beta_term[d] = fh->beta*(arr[C].u_fh[d]*arr[C].ro_fh - arr[C].uro_md[d]);
                }
                else if(fh->scheme == Two_Way)
                {
//                  arr[C].delta_ro = arr[C].ron_prime;
//                  for(int d = 0; d < DIM; d++)
//                      arr[C].beta_term[d] = fh->beta*arr[C].mn_prime[d];
                    arr[C].delta_ro = arr[C].ro_fh - arr[C].ro_md;                                          // Layer n may work better than n+1/2
                    for(int d = 0; d < DIM; d++)
                    {
                        arr[C].beta_term[d] = fh->beta*(arr[C].u_fh[d]*arr[C].ro_fh - arr[C].uro_md[d]);    // Layer n may work better than n+1/2

                        //arr[C].alpha_term_exp[d] = fh->alpha*(arr[C].S*(1 - arr[C].S) + arr[C].S)*arr[C].delta_ro*arr[C].inv_ro*fh->grid.h[C][d];

                        switch(d)
                        {
                        case 0:
                            ij = i;
                            break;
                        case 1:
                            ij = j;
                            break;
                        case 2:
                            ij = k;
                            break;
                        }
                        ij -= fh->N_shift[d];
                        switch(ij)
                        {
                        case 0:
                            arr[C].alpha_term_exp[d]  = -fh->alpha*arr[C].delta_ro*arr[C].inv_ro*fh->grid.h[C][d]*0.5;
                            arr[C].alpha_term_exp[d] += -fh->alpha*(arr[CR].ro_fh - arr[CR].ro_md)/arr[CR].ro_md*fh->grid.h[CR][d];
                            break;
                        case 1:
                            arr[C].alpha_term_exp[d] = -fh->alpha*arr[C].delta_ro*arr[C].inv_ro*fh->grid.h[C][d]*0.5;
                            break;
                        case 2:
                            arr[C].alpha_term_exp[d] = 0;
                            break;
                        case 3:
                            arr[C].alpha_term_exp[d] = fh->alpha*arr[C].delta_ro*arr[C].inv_ro*fh->grid.h[C][d]*0.5;
                            break;
                        case 4:
                            arr[C].alpha_term_exp[d]  = fh->alpha*arr[C].delta_ro*arr[C].inv_ro*fh->grid.h[C][d]*0.5;
                            arr[C].alpha_term_exp[d] += fh->alpha*(arr[CL].ro_fh - arr[CL].ro_md)/arr[CL].ro_md*fh->grid.h[CL][d];
                            break;
                        }
                    }

                    if(fh->grid.md[C] == FH_zone) arr[C].delta_ro = 0;
                }
            }
        }
    }

    for(int k = fh->N_shift[2]; k < (fh->N_md[2] + fh->N_shift[2]); k++)
    {
        for(int j = fh->N_shift[1]; j < (fh->N_md[1] + fh->N_shift[1]); j++)
        {
            for(int i = fh->N_shift[0]; i < (fh->N_md[0] + fh->N_shift[0]); i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                {
                    arr[Cm].grad_ro[d] = fh->alpha*(arr[CRm].delta_ro - arr[CLm].delta_ro)/(0.5*(fh->grid.h[CLm][d] + 2.0*fh->grid.h[Cm][d] + fh->grid.h[CRm][d]));

                    for(int du = 0; du < DIM; du++)
                        //arr[Cm].alpha_u_grad[du][d] = arr[Cm].grad_ro[d]*arr[Cm].S*(1 - arr[Cm].S)*arr[Cm].u_md[du];    // TODO: Fast but rough estimation!
                        arr[Cm].alpha_u_grad[du][d] = arr[Cm].alpha_term_exp[d]*arr[Cm].uro_md[du];
                }
            }
        }
    }

    for(int k = fh->N_shift[2]; k < (fh->N_md[2] + fh->N_shift[2]); k++)
    {
        for(int j = fh->N_shift[1]; j < (fh->N_md[1] + fh->N_shift[1]); j++)
        {
            for(int i = fh->N_shift[0]; i < (fh->N_md[0] + fh->N_shift[0]); i++)
            {
                ASSIGN_IND(ind, i, j, k);

                if((fh->scheme == Two_Way) && (fh->grid.md[C] == boundary)) continue;

                for(int du = 0; du < DIM; du++)
                {
                    for(int d = 0; d < DIM; d++)
                    {
                        alpha_term[d] = (arr[CRm].alpha_u_grad[du][d] - arr[CLm].alpha_u_grad[du][d])
                                /(0.5*(fh->grid.h[CLm][d] + 2.0*fh->grid.h[Cm][d] + fh->grid.h[CRm][d]));
                    }

                    arr[Cm].alpha_term[du] = SUM(alpha_term);
                }
            }
        }
    }
}

