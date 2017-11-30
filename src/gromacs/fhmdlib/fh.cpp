#include "data_structures.h"
#include "macro.h"


void define_FH_grid(t_commrec *cr, FHMD *fh)
{
    dvec h0;
    ivec ind;

    for(int d = 0; d < DIM; d++)
        h0[d] = fh->box[d]/(double)(fh->N[d]);

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                {
                    fh->grid.h[C][d] = h0[d];
                    fh->grid.n[C][d] = (double)(ind[d])*h0[d];
                    fh->grid.c[C][d] = fh->grid.n[C][d] + 0.5*h0[d];
                }
            }
        }
    }

    for(int i = 0; i < fh->Ntot; i++)
    {
        fh->grid.vol[i]  = fh->grid.h[i][0]*fh->grid.h[i][1]*fh->grid.h[i][2];
        fh->grid.ivol[i] = 1.0/fh->grid.vol[i];
    }

#ifdef FHMD_GRID_DEBUG
    if(MASTER(cr)) {
        printf(MAKE_YELLOW "FHMD DEBUG: Grid X (i, nodes, centres, steps):\n");
        for(int i = 0; i < fh->N[0]; i++)
            printf("%3d  %8.4f  %8.4f  %8.4f\n", i, fh->grid.n[I3(i,0,0,fh->N)][0], fh->grid.c[I3(i,0,0,fh->N)][0], fh->grid.h[I3(i,0,0,fh->N)][0]);
        printf("%3d  %8.4f     ----      ----\n", fh->N[0], fh->grid.n[I3(fh->N[0],0,0,fh->N)][0]);
        printf(MAKE_YELLOW "FHMD DEBUG: Grid Y (i, nodes, centres, steps):\n");
        for(int i = 0; i < fh->N[1]; i++)
            printf("%3d  %8.4f  %8.4f  %8.4f\n", i, fh->grid.n[I3(0,i,0,fh->N)][1], fh->grid.c[I3(0,i,0,fh->N)][1], fh->grid.h[I3(0,i,0,fh->N)][1]);
        printf("%3d  %8.4f     ----      ----\n", fh->N[1], fh->grid.n[I3(0,fh->N[1],0,fh->N)][1]);
        printf(MAKE_YELLOW "FHMD DEBUG: Grid Z (i, nodes, centres, steps):\n");
        for(int i = 0; i < fh->N[2]; i++)
            printf("%3d  %8.4f  %8.4f  %8.4f\n", i, fh->grid.n[I3(0,0,i,fh->N)][2], fh->grid.c[I3(0,0,i,fh->N)][2], fh->grid.h[I3(0,0,i,fh->N)][2]);
        printf("%3d  %8.4f     ----      ----\n", fh->N[2], fh->grid.n[I3(0,0,fh->N[2],fh->N)][2]);
        printf("==============================================\n");
        printf(RESET_COLOR "\n");
    }
#endif
}
