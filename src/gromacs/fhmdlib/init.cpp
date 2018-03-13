#include "data_structures.h"
#include "parser.h"
#include "fh_functions.h"
#include "estimate.h"
#include "sfunction.h"


int fhmd_init(matrix box, int N_atoms, real mass[], rvec x[], double dt_md, gmx_mtop_t *mtop, t_commrec *cr, FHMD *fh)
{
    int N_atoms_th = N_atoms;

    FILE *fw;

    if(MASTER(cr))
    {
        char const *fname_in  = "coupling.prm";
        char const *fname_out = "coupling_out.prm";

        /* Initial output */

        printf(MAKE_GREEN "\n  Aston University, Department of Mathematics, Dmitry Nerukh Research Group\n");
        printf("  Queen Mary University of London, School of Engineering and Material Science\n\n");
        printf(MAKE_PURPLE "     Hybrid Molecular Dynamics - 2-Way Coupling Parallel VERSION %4.2f\n\n", FHMD_VERSION);

        /* Default values of FHMD parameters */

        fh->scheme      = 1;
        fh->S           = 0;
        fh->R1          = 0.5;
        fh->R2          = 0.9;
        fh->Smin        = 0;
        fh->Smax        = 0.99;
        fh->alpha       = 100;
        fh->beta        = 100;
        fh->gamma_x     = 0;
        fh->gamma_u     = 0;
        fh->eps_rho     = 0;
        fh->eps_mom     = 0;
        fh->S_berendsen = 1;
        fh->N[0]        = 5;
        fh->N[1]        = 5;
        fh->N[2]        = 5;
        fh->N_md[0]     = fh->N[0];
        fh->N_md[1]     = fh->N[1];
        fh->N_md[2]     = fh->N[2];
        fh->FH_EOS      = 1;
        fh->FH_step     = 10;
        fh->FH_equil    = 10000;
        fh->FH_dens     = 600;
        fh->FH_temp     = 298.15;
        fh->FH_blend    = 0.005;
        fh->Noutput     = 10;

        /* Read FHMD parameters */

        printf(MAKE_GREEN "FHMD: Reading parameters from %s..." RESET_COLOR " ", fname_in);

        int ok = parse_prm(fname_in, fh);

        if(ok == 1) {
            printf(MAKE_GREEN "...OK\n" RESET_COLOR "\n");
            fw = fopen(fname_out, "w");
        } else if(ok == -1) {
            printf(MAKE_RED "\nFHMD: File %s not found. Generating default parameter file...\n" RESET_COLOR "\n", fname_in);
            fw = fopen(fname_in, "w");
        } else {
            printf(MAKE_RED "\nFHMD: ERROR in %s file\n" RESET_COLOR "\n", fname_in);
            exit(2);
        }

        /* Print FHMD parameters to the screen and output file */

        fprintf(fw, "; Hybrid Molecular Dynamics - 2-Way Coupling Parallel VERSION %4.2f\n\n", FHMD_VERSION);

        fprintf(fw, "scheme = %d              ; 0 - Pure MD, 1 - One-way coupling, 2 - Two-way coupling\n\n", fh->scheme);

        fprintf(fw, "S = %g                   ; Parameter S (-1 - fixed sphere, -2 - moving sphere)\n\n", fh->S);
        fprintf(fw, "R1   = %g              ; MD sphere radius for variable S, [0..1]\n", fh->R1);
        fprintf(fw, "R2   = %g              ; FH sphere radius for variable S, [0..1]\n", fh->R2);
        fprintf(fw, "Smin = %g                ; Minimum S for variable S\n", fh->Smin);
        fprintf(fw, "Smax = %g             ; Maximum S for variable S\n\n", fh->Smax);

        switch(fh->scheme)
        {
        case Pure_MD:
            printf(MAKE_RED "FHMD: Starting Pure MD simulation (without MD/FH coupling)\n" RESET_COLOR "\n");
            break;
        case One_Way:
            printf(MAKE_YELLOW "FHMD: One-way MD/FH coupling\n" RESET_COLOR "\n");
            break;
        case Two_Way:
            break;
        }

        if(fh->S >= 0.0)
        {
            printf(MAKE_PURPLE "FHMD: S = %g\n", fh->S);
            fh->S_function = constant_S;
        }
        else
        {
            printf(MAKE_PURPLE "FHMD: S = S(x,y,z) = [%g, %g] with R1 = %g, R2 = %g\n", fh->Smin, fh->Smax, fh->R1, fh->R2);
            fh->S_function = fixed_sphere;

            fh->R1 *= box[0][0]*0.5;
            fh->R2 *= box[0][0]*0.5;
            fh->R12 = fh->R1*fh->R1;
            fh->R22 = fh->R2*fh->R2;
            fh->RS  = (fh->Smax - fh->Smin)/(fh->R2 - fh->R1);
            printf(MAKE_GREEN "FHMD: Absolute values of R [nm]: R1 = %f, R2 = %f\n", fh->R1, fh->R2);

            if(fh->S < -1.5)
            {
                printf(MAKE_PURPLE "FHMD: The MD/FH sphere will follow the protein\n");
                fh->S_function = moving_sphere;
            }
        }

        printf(MAKE_GREEN "FHMD: alpha = %g [nm^2/ps], beta = %g [ps^-1]\n", fh->alpha, fh->beta);
        fprintf(fw, "alpha   = %g           ; Alpha parameter for dx/dt and du/dt equations, nm^2/ps\n", fh->alpha);
        fprintf(fw, "beta    = %g           ; Beta parameter for du/dt equation, ps^-1\n\n", fh->beta);

        if(fh->scheme == Two_Way)
        {
            printf("FHMD: MD dissipator: gamma_x = %g [ps^-1], gamma_u = %g [ps^-1]\n", fh->gamma_x, fh->gamma_u);
            printf("FHMD: FH dissipator: eps_rho = %g [--], eps_mom = %g [--]\n", fh->eps_rho, fh->eps_mom);
            fprintf(fw, "gamma_x = %g             ; Gamma_x parameter (MD density fluctuations dissipator), ps^-1\n", fh->gamma_x);
            fprintf(fw, "gamma_u = %g             ; Gamma_u parameter (MD velocity fluctuations dissipator), ps^-1\n", fh->gamma_u);
            fprintf(fw, "eps_rho = %g             ; Eps_rho parameter (FH density fluctuations dissipator)\n", fh->eps_rho);
            fprintf(fw, "eps_mom = %g             ; Eps_mom parameter (FH momentum fluctuations dissipator)\n\n", fh->eps_mom);
        }

        if(fh->S_berendsen >= 0)
            printf("FHMD: Berendsen thermostat works for S <= %g\n", fh->S_berendsen);
        else
            printf("FHMD: Berendsen thermostat with (1 - S^%g) multiplier\n", -fh->S_berendsen);
        fprintf(fw, "S_berendsen = %g         ; If S_berendsen >= 0, Berendsen thermostat will work for S <= S_berendsen,\n", fh->S_berendsen);
        fprintf(fw, "                        ; otherwise factor (1-S^(-S_berendsen)) will be applied (local thermostat)\n\n");

        for(int d = 0; d < DIM; d++)
        {
            fh->box[d]         = box[d][d];
            fh->box05[d]       = 0.5*fh->box[d];
            fh->protein_com[d] = fh->box05[d];
            fh->N_shift[d]     = (fh->N[d] - fh->N_md[d])/2;
        }

        fh->box_volume = fh->box[0]*fh->box[1]*fh->box[2];

        printf("FHMD: MD box size:  %g x %g x %g [nm]\n", fh->box[0], fh->box[1], fh->box[2]);
        printf("FHMD: FH grid size: %d x %d x %d\n", fh->N[0], fh->N[1], fh->N[2]);
        fprintf(fw, "Nx = %d                  ; Number of FH cells along X axis\n", fh->N[0]);
        fprintf(fw, "Ny = %d                  ; Number of FH cells along Y axis\n", fh->N[1]);
        fprintf(fw, "Nz = %d                  ; Number of FH cells along Z axis\n\n", fh->N[2]);

        if(fh->scheme == Two_Way)
        {
            printf("FHMD: Small-scale MD/FH grid size: %d x %d x %d\n", fh->N_md[0], fh->N_md[1], fh->N_md[2]);
            fprintf(fw, "NxMD = %d                ; Number of small-scale MD-FH cells along X axis\n", fh->N_md[0]);
            fprintf(fw, "NyMD = %d                ; Number of small-scale MD-FH cells along Y axis\n", fh->N_md[1]);
            fprintf(fw, "NzMD = %d                ; Number of small-scale MD-FH cells along Z axis\n\n", fh->N_md[2]);
        }
        else
        {
            for(int d = 0; d < DIM; d++)
            {
                fh->N_shift[d] = 0;
                fh->N_md[d]    = fh->N[d];
            }
        }

        fh->Ntot    = fh->N[0]*fh->N[1]*fh->N[2];
        fh->Ntot_md = fh->N_md[0]*fh->N_md[1]*fh->N_md[2];

        switch(fh->FH_EOS)
        {
        case 0:
            fh->eos = eos_argon;
            printf("FHMD: Equation of state: Liquid Argon (300K)\n");
            break;
        case 1:
            fh->eos = eos_spce;
            printf("FHMD: Equation of state: Rigid SPC/E water\n");
            break;
        default:
            printf(MAKE_RED "\nFHMD: Unknown equation of state (%d) in %s\n" RESET_COLOR "\n", fh->FH_EOS, fname_in);
            exit(18);
        }

        fprintf(fw, "FH_EOS   = %d            ; EOS: 0 - Liquid Argon, 1 - SPC/E water\n", fh->FH_EOS);

        fh->dt_FH = (double)(fh->FH_step)*dt_md;

        printf("FHMD: FH time step dt_FH = %d * dt_MD = %g [ps]\n", fh->FH_step, fh->dt_FH);
        fprintf(fw, "FH_step  = %d           ; FH time step dt_FH = FH_step * dt_MD\n", fh->FH_step);

        if(fh->scheme == Two_Way) fh->FH_equil = 0;
        printf("FHMD: FH equilibration steps: %d\n", fh->FH_equil);
        fprintf(fw, "FH_equil = %d        ; Number of time steps for the FH model equilibration (for 1-way coupling)\n", fh->FH_equil);

        printf("FHMD: FH Density = %g [amu/nm^3], FH Temperature = %g [K]\n", fh->FH_dens, fh->FH_temp);
        fprintf(fw, "FH_dens  = %g          ; FH mean density\n", fh->FH_dens);
        fprintf(fw, "FH_temp  = %g       ; FH mean temperature\n", fh->FH_temp);
        fprintf(fw, "FH_blend = %g        ; FH Blending: -1 - dynamic, or define static blending parameter (0..1)\n\n", fh->FH_blend);

        printf("FHMD: MD/FH arrays will be written every %d MD time steps\n", fh->Noutput);
        fprintf(fw, "Noutput  = %d           ; Write arrays to files every Noutput MD time steps (0 - do not write)\n", fh->Noutput);

        printf(RESET_COLOR "\n");

        fflush(stdout);
    } // if(MASTER(cr))

    fhmd_reset_statistics(fh);

    fh->total_density = 0;
    for(int i = 0; i < N_atoms_th; i++)
        fh->total_density += mass[i];

    /* Broadcast parameters to all threads */

    if(PAR(cr))
    {
        gmx_sumd(1, &fh->total_density, cr);
        gmx_sumi(1, &N_atoms, cr);
        gmx_bcast(sizeof(FHMD), fh, cr);
    }

    fh->total_density /= fh->box_volume;

    if(fh->eos == eos_spce)
    {
        fhmd_find_protein(mtop, N_atoms_th, mass, cr, fh);
        if(fh->protein_N)
            fhmd_find_protein_com(mtop, N_atoms_th, x, mass, cr, fh);
    }
    else
    {
        fh->protein_N    = 0;
        fh->protein_mass = 0;
    }

    if(MASTER(cr))
    {
        printf(MAKE_GREEN "FHMD: Total number of atoms in the box: %d\n", N_atoms);
        printf("FHMD: Total density of the box: %g [amu/nm^3]\n", fh->total_density);

        if(fh->protein_N > 0)
            printf(MAKE_PURPLE "FHMD: Found protein: %d atoms, mass = %g [amu], COM = (%g, %g, %g) [nm]\n",
                    fh->protein_N, fh->protein_mass, fh->protein_com[0], fh->protein_com[1], fh->protein_com[2]);

        printf(RESET_COLOR "\n");
        fflush(stdout);

        fprintf(fw, "\n; You may consider to use FH_dens = %g since this is the total MD density of the box.\n", fh->total_density);
        fprintf(fw, "; NB: Please use spaces before and after '=' in this file, e.g. 'S = 0' (not 'S=0').\n");
        fclose(fw);
    }

    if(fh->scheme == Pure_MD) return 0;     // Start Pure MD simulation

    /* Allocate memory */

    fh->arr  = (FH_arrays*)calloc(fh->Ntot, sizeof(FH_arrays));
    fh->ind  = (int*)calloc(N_atoms, sizeof(int));
    fh->indv = (ivec*)calloc(N_atoms, sizeof(ivec));
    fh->vel  = (dvec*)calloc(N_atoms, sizeof(dvec));

    fh->mpi_linear = (double*)malloc(8*fh->Ntot*sizeof(double));   // 8 components: ro_md, uro_md[3], ro_md_s, uro_md_s[3]

    if(fh->arr == NULL || fh->ind == NULL || fh->indv == NULL || fh->vel == NULL || fh->mpi_linear == NULL)
    {
        if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (array allocator)\n" RESET_COLOR "\n");
        fflush(stdout);
        exit(3);
    }

    fh->grid.c    = (dvec*)malloc(fh->Ntot*sizeof(dvec));
    fh->grid.n    = (dvec*)malloc(fh->Ntot*sizeof(dvec));
    fh->grid.h    = (dvec*)malloc(fh->Ntot*sizeof(dvec));
    fh->grid.vol  = (double*)malloc(fh->Ntot*sizeof(double));
    fh->grid.ivol = (double*)malloc(fh->Ntot*sizeof(double));
    fh->grid.md   = (FHMD_CELL*)malloc(fh->Ntot*sizeof(FHMD_CELL));

    if(fh->grid.c == NULL || fh->grid.n == NULL || fh->grid.h == NULL || fh->grid.vol == NULL || fh->grid.ivol == NULL)
    {
        if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (FH grid allocator)\n" RESET_COLOR "\n");
        fflush(stdout);
        exit(3);
    }

    fh->stat.avg_rho_md_cell = (double*)calloc(fh->Ntot, sizeof(double));
    fh->stat.avg_rho_fh_cell = (double*)calloc(fh->Ntot, sizeof(double));

    if(fh->stat.avg_rho_md_cell == NULL || fh->stat.avg_rho_fh_cell == NULL)
    {
        if(MASTER(cr)) printf(MAKE_RED "\nFHMD: ERROR: Out of memory (Statistics allocator)\n" RESET_COLOR "\n");
        fflush(stdout);
        exit(3);
    }

    /* Create FH grid and initialize FH solver */

    define_FH_grid(cr, fh);
    if(fh->scheme == One_Way)
        FH_init(fh, cr);

    for(int i = 0; i < fh->Ntot; i++)
    {
        fh->arr[i].S = 1;
        for(int d = 0; d < DIM; d++)
        {
            fh->arr[i].Sf[d] = 1;
        }
    }

    for(int i = 0; i < FHMD_COUETTE_LAYERS; i++)
    {
        fh->avg_vel_tot[i]   = 0;
        fh->avg_vel_S_tot[i] = 0;
        fh->avg_n_tot[i]     = 0;
    }

    if(MASTER(cr))
    {
        if(fh->scheme == One_Way)
        {
            FH_equilibrate(fh);
            printf(MAKE_GREEN "FHMD: Initialization finished. Starting MD/FH solver...\n" RESET_COLOR "\n");
        }
        fflush(stdout);
    }

    return 1;   // Success
}
