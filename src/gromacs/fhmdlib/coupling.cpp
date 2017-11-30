#include "data_structures.h"
#include "macro.h"


void fhmd_update_MD_in_FH(rvec x[], rvec v[], real mass[], int N_atoms, FHMD *fh)
{

    FH_arrays *arr = fh->arr;
    dvec       xn;
    int        ind;

    /* Reset statistics */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md = 0;
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

        if(ind < 0)     // This should never happen... only if the coordinates are NaN.
        {
            printf(MAKE_RED "\nFHMD: ERROR: Solution diverged. Atom's coordinates: (%g, %g, %g) nm\n" RESET_COLOR "\n", xn[0], xn[1], xn[2]);
            exit(21);
        }

        fh->ind[n] = ind;

        arr[ind].ro_md += mass[n];
    }

    /* Update statistics */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md *= fh->grid.ivol[i];
    }

}


void fhmd_sum_arrays(t_commrec *cr, FHMD *fh)
{

    FH_arrays *arr = fh->arr;

    /* Pack FHMD arrays to linear array */
    for(int i = 0; i < fh->Ntot; i++)
    {
        fh->mpi_linear[i] = arr[i].ro_md;
    }

    /* Broadcast linear array */
    gmx_sumd(fh->Ntot, fh->mpi_linear, cr);

    /* Unpack linear array */
    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md = fh->mpi_linear[i];
    }

}


void fhmd_calculate_MDFH_terms(FHMD *fh)
{

    FH_arrays *arr = fh->arr;

    for(int i = 0; i < fh->Ntot; i++)
    {
        if(arr[i].ro_md <= 0) {
            printf(MAKE_RED "\nFHMD: ERROR: Zero or NaN MD density in the cell #%d (ro_md = %g)\n" RESET_COLOR "\n", i, arr[i].ro_md);
            exit(22);
        }
        arr[i].inv_ro = 1.0/arr[i].ro_md;
    }

}

