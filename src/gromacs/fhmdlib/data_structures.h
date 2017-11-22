#ifndef FHMD_DATA_STRUCTURES_H_
#define FHMD_DATA_STRUCTURES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "params.h"
#include "gromacs/mdtypes/commrec.h"        /* GROMACS MPI definitions, 't_commrec', MASTER(), PAR(), etc. */
#include "gromacs/gmxlib/network.h"         /* GROMACS MPI functions, gmx_bcast(), gmx_sumf(), etc. */


typedef struct FH_VecI          /* Integer vector */
{
    int x, y, z;
} FH_VecI;


typedef struct FH_VecD          /* Double vector */
{
    double x, y, z;
} FH_VecD;


typedef struct FH_arrays        /* FH/MD arrays */
{
    double      ro_md, ro_fh;   /* densities */
} FH_arrays;


typedef struct FH_grid          /* Computational grid */
{
    FH_VecD    *c;              /* FH cell centres coordinates */
    FH_VecD    *n;              /* FH cell nodes coordinates */
    FH_VecD    *h;              /* FH cell steps */
    double     *vol;            /* FH cell volume */
    double     *ivol;           /* 1/cellVolume */
} FH_grid;


typedef struct FHMD
{

    FH_arrays  *arr;            /* FH/MD arrays */
    FH_grid     grid;           /* FH grid */
    int        *ind;            /* FH cell number for each atom */
    double     *mpi_linear;     /* Linear array to synchronize MDFH arrays */

    double      S;              /* Parameter S (-1 if S is variable) */
    double      R1;             /* MD sphere radius for variable S, [0..1] */
    double      R2;             /* FH sphere radius for variable S, [0..1] */
    double      Smin;           /* Minimum S for variable S */
    double      Smax;           /* Maximum S for variable S */
    double      R12, R22, RS;   /* Derived variables from R1, R2, Smin, Smax */
    double      alpha;          /* Alpha parameter for dx/dt and du/dt equations, nm^2/ps */
    double      beta;           /* Beta parameter, nm^2/ps or ps^-1 depending on the scheme */

    FH_VecI     N;              /* Number of FH cells along each direction */
    FH_VecD     box;            /* Box size */
    int         Ntot;           /* Total number of FH cells */

    int         FH_step;        /* dt_FH = FH_step * dt_MD */
    int         FH_equil;       /* Number of time steps for the FH model equilibration */
    double      FH_dens;        /* FH mean density */
    double      FH_temp;        /* FH mean temperature */

} FHMD;

#endif /* FHMD_DATA_STRUCTURES_H_ */
