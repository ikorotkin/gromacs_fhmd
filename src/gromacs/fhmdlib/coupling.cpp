#include "data_structures.h"
#include "macro.h"


void fhmd_update_MD_in_FH(rvec x[], rvec v[], real mass[], int N_atoms, FHMD *fh)
{

    FH_arrays *arr = fh->arr;
    FH_VecI    ind;
    int        index;

    /* Reset statistics */

    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_md = 0;
    }

    /* Collect statistics */

    for(int n = 0; n < N_atoms; n++)
    {

        ind.x = (int)(x[n][0]/fh->box.x*(double)(fh->N.x));
        ind.y = (int)(x[n][1]/fh->box.y*(double)(fh->N.y));
        ind.z = (int)(x[n][2]/fh->box.z*(double)(fh->N.z));

        index = I(ind.x, ind.y, ind.z, fh->N);

        fh->ind[n] = index;

        arr[index].ro_md += mass[n];

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
