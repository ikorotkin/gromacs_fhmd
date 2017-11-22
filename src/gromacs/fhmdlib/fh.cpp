#include "data_structures.h"
#include "macro.h"

#define FHMD_GRID_DEBUG


void define_FH_grid(t_commrec *cr, FHMD *fh)
{
    FH_VecD h0;

    h0.x = (double)(fh->box.x)/(double)(fh->N.x);
    h0.y = (double)(fh->box.y)/(double)(fh->N.y);
    h0.z = (double)(fh->box.z)/(double)(fh->N.z);

    /* Uniform grid */

    for(int k = 0; k < fh->N.z; k++)
    {
        for(int j = 0; j < fh->N.y; j++)
        {
            for(int i = 0; i < fh->N.x; i++)
            {
                fh->grid.h[C].x = h0.x;
                fh->grid.h[C].y = h0.y;
                fh->grid.h[C].z = h0.z;

                fh->grid.n[C].x = (double)(i)*h0.x;
                fh->grid.n[C].y = (double)(j)*h0.y;
                fh->grid.n[C].z = (double)(k)*h0.z;

                fh->grid.c[C].x = fh->grid.n[C].x + 0.5*h0.x;
                fh->grid.c[C].y = fh->grid.n[C].y + 0.5*h0.y;
                fh->grid.c[C].z = fh->grid.n[C].z + 0.5*h0.z;
            }
        }
    }

    for(int i = 0; i < fh->Ntot; i++)
    {
        fh->grid.vol[i]  = fh->grid.h[i].x*fh->grid.h[i].y*fh->grid.h[i].z;
        fh->grid.ivol[i] = 1.0/fh->grid.vol[i];
    }

#ifdef FHMD_GRID_DEBUG
    if(MASTER(cr)) {
        printf(MAKE_YELLOW "FHMD DEBUG: Grid X (i, nodes, centres, steps):\n");
        for(int i = 0; i < fh->N.x; i++)
            printf("%3d  %8.4f  %8.4f  %8.4f\n", i, fh->grid.n[I(i,0,0,fh->N)].x, fh->grid.c[I(i,0,0,fh->N)].x, fh->grid.h[I(i,0,0,fh->N)].x);
        printf("%3d  %8.4f     ----      ----\n", fh->N.x, fh->grid.n[I(fh->N.x,0,0,fh->N)].x);
        printf(MAKE_YELLOW "FHMD DEBUG: Grid Y (i, nodes, centres, steps):\n");
        for(int i = 0; i < fh->N.y; i++)
            printf("%3d  %8.4f  %8.4f  %8.4f\n", i, fh->grid.n[I(0,i,0,fh->N)].y, fh->grid.c[I(0,i,0,fh->N)].y, fh->grid.h[I(0,i,0,fh->N)].y);
        printf("%3d  %8.4f     ----      ----\n", fh->N.x, fh->grid.n[I(0,fh->N.y,0,fh->N)].y);
        printf(MAKE_YELLOW "FHMD DEBUG: Grid Z (i, nodes, centres, steps):\n");
        for(int i = 0; i < fh->N.z; i++)
            printf("%3d  %8.4f  %8.4f  %8.4f\n", i, fh->grid.n[I(0,0,i,fh->N)].z, fh->grid.c[I(0,0,i,fh->N)].z, fh->grid.h[I(0,0,i,fh->N)].z);
        printf("%3d  %8.4f     ----      ----\n", fh->N.x, fh->grid.n[I(0,0,fh->N.z,fh->N)].z);
        printf("==============================================\n");
        printf(RESET_COLOR "\n");
    }
#endif

}
