#include "data_structures.h"
#include "parser.h"
#include "fh.h"


int fhmd_init(matrix box, int N_atoms, t_commrec *cr, FHMD *fh)
{

    if(MASTER(cr))
    {

        char const *fname_in  = "coupling.prm";
        char const *fname_out = "coupling_out.prm";

        /* Initial output */

        printf(MAKE_GREEN "\n  Aston University, Department of Mathematics, Dmitry Nerukh Research Group\n");
        printf("  Queen Mary University of London, School of Engineering and Material Science\n\n");
        printf(MAKE_BLUE "     Hybrid Molecular Dynamics - 2-Way Coupling Parallel VERSION %4.2f\n\n", fhmd_version);

        /* Default FHMD parameters */

        fh->S           = 0;
        fh->R1          = 0.5;
        fh->R2          = 0.9;
        fh->Smin        = 0;
        fh->Smax        = 0.99;
        fh->alpha       = 10;
        fh->beta        = 10;
        fh->N.x         = 5;
        fh->N.y         = 5;
        fh->N.z         = 5;
        fh->FH_step     = 10;
        fh->FH_equil    = 1000;
        fh->FH_dens     = 600;
        fh->FH_temp     = 300;

        /* Read FHMD parameters */

        printf(MAKE_GREEN "FHMD: Reading parameters from %s..." RESET_COLOR " ", fname_in);

        int ok = parse_prm(fname_in, fh);
        
        FILE *fw;

        if(ok == 1) {
            printf(MAKE_GREEN "...OK\n" RESET_COLOR "\n");
            fw = fopen(fname_out, "w");
        } else if(ok == -1) {
            printf(MAKE_RED "\nFHMD: File %s not found. Generating one with default parameters...\n" RESET_COLOR "\n", fname_in);
            fw = fopen(fname_in, "w");
        } else {
            printf(MAKE_RED "\nFHMD: ERROR in %s file\n" RESET_COLOR "\n", fname_in);
            exit(2);
        }

        /* Print FHMD parameters to the screen and output file */

        fprintf(fw, "; Hybrid Molecular Dynamics - 2-Way Coupling Parallel VERSION %4.2f\n\n", fhmd_version);

        fprintf(fw, "S = %g               ; Parameter S (-1 if S is variable)\n\n", fh->S);
        fprintf(fw, "R1   = %g          ; MD sphere radius for variable S, [0..1]\n", fh->R1);
        fprintf(fw, "R2   = %g          ; FH sphere radius for variable S, [0..1]\n", fh->R2);
        fprintf(fw, "Smin = %g            ; Minimum S for variable S\n", fh->Smin);
        fprintf(fw, "Smax = %g         ; Maximum S for variable S\n\n", fh->Smax);

        if(fh->S >= 0.0) {
            printf(MAKE_PURPLE "FHMD: S = %g\n", fh->S);
        } else {
            printf(MAKE_PURPLE "FHMD: S = S(x,y,z) = [%g, %g] with R1 = %g, R2 = %g\n", fh->Smin, fh->Smax, fh->R1, fh->R2);
            fh->R1 *= box[0][0]*0.5;
            fh->R2 *= box[0][0]*0.5;
            fh->R12 = fh->R1*fh->R1;
            fh->R22 = fh->R2*fh->R2;
            fh->RS  = (fh->Smax - fh->Smin)/(fh->R2 - fh->R1);
            printf(MAKE_GREEN "FHMD: Absolute values of R [nm]: R1 = %f, R2 = %f\n", fh->R1, fh->R2);
        }

        printf(MAKE_GREEN "FHMD: alpha = %g [nm^2/ps], beta = %g [nm^2/ps or ps^-1]\n", fh->alpha, fh->beta);
        fprintf(fw, "alpha = %g          ; Alpha parameter for dx/dt and du/dt equations, nm^2/ps\n", fh->alpha);
        fprintf(fw, "beta  = %g          ; Beta parameter, nm^2/ps or ps^-1 depending on the scheme\n\n", fh->beta);

        fh->box.x = box[0][0];
        fh->box.y = box[1][1];
        fh->box.z = box[2][2];

        printf("FHMD: MD/FH box size: %g x %g x %g\n", fh->box.x, fh->box.y, fh->box.z);
        printf("FHMD: FH grid size:   %d x %d x %d\n", fh->N.x, fh->N.y, fh->N.z);
        fprintf(fw, "Nx = %d              ; Number of FH cells along X axis\n", fh->N.x);
        fprintf(fw, "Ny = %d              ; Number of FH cells along Y axis\n", fh->N.y);
        fprintf(fw, "Nz = %d              ; Number of FH cells along Z axis\n\n", fh->N.z);

        fh->Ntot  = fh->N.x*fh->N.y*fh->N.z;

        printf("FHMD: FH time step dt_FH = %d * dt_MD\n", fh->FH_step);
        fprintf(fw, "FH_step  = %d       ; FH time step dt_FH = FH_step * dt_MD\n", fh->FH_step);

        printf("FHMD: FH equilibration steps: %d\n", fh->FH_equil);
        fprintf(fw, "FH_equil = %d     ; Number of time steps for the FH model equilibration\n", fh->FH_equil);

        printf("FHMD: FH Density = %g [amu/nm^3], FH Temperature = %g [K]\n", fh->FH_dens, fh->FH_temp);
        fprintf(fw, "FH_dens  = %g      ; FH mean density\n", fh->FH_dens);
        fprintf(fw, "FH_temp  = %g      ; FH mean temperature\n", fh->FH_temp);

        printf(RESET_COLOR "\n");

        fflush(stdout);
        fclose(fw);

        /* Open files for writing */


    } // if(MASTER(cr))

    /* Broadcast parameters to all threads */

    if(PAR(cr))
    {
        gmx_bcast(sizeof(FHMD), fh, cr);
        gmx_sumi(1, &N_atoms, cr);
    }

    if(MASTER(cr))
    {
        printf(MAKE_GREEN "FHMD: Total number of atoms in the box: %d\n", N_atoms);
        printf(RESET_COLOR "\n");
        fflush(stdout);
    }

    /* Allocate memory */

    fh->arr = (FH_arrays*)calloc(fh->Ntot, sizeof(FH_arrays));
    fh->ind = (int*)calloc(N_atoms, sizeof(int));

    fh->mpi_linear = (double*)malloc(fh->Ntot*sizeof(double));

    if(fh->arr == NULL || fh->ind == NULL || fh->mpi_linear == NULL)
    {
        if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (array allocator)\n" RESET_COLOR "\n");
        fflush(stdout);
        exit(3);
    }

    fh->grid.c    = (FH_VecD*)malloc(fh->Ntot*sizeof(FH_VecD));
    fh->grid.n    = (FH_VecD*)malloc(fh->Ntot*sizeof(FH_VecD));
    fh->grid.h    = (FH_VecD*)malloc(fh->Ntot*sizeof(FH_VecD));
    fh->grid.vol  = (double*)malloc(fh->Ntot*sizeof(double));
    fh->grid.ivol = (double*)malloc(fh->Ntot*sizeof(double));

    if(fh->grid.c == NULL || fh->grid.n == NULL || fh->grid.h == NULL || fh->grid.vol == NULL || fh->grid.ivol == NULL)
    {
        if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (FH grid allocator)\n" RESET_COLOR "\n");
        fflush(stdout);
        exit(3);
    }

    /* Create FH grid */

    define_FH_grid(cr, fh);

    return 1;   // Success

}
