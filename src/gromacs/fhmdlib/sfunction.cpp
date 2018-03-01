#include "data_structures.h"
#include "macro.h"

#include "gromacs/topology/mtop_util.h"     /* This is for gmx_mtop_atominfo_global() */


double fhmd_Sxyz_r(const rvec x, const dvec c, FHMD *fh)
{
    dvec xd;

    for(int d = 0; d < DIM; d++)
    {
        xd[d] = fabs(x[d] - c[d]);
        if(xd[d] > fh->box05[d]) xd[d] -= fh->box[d];
    }

    double r2 = xd[0]*xd[0] + xd[1]*xd[1] + xd[2]*xd[2];

    if(r2 <= fh->R12) return fh->Smin;      // S = Smin inside the radius R1
    if(r2 >= fh->R22) return fh->Smax;      // S = Smax outside the radius R2

    double r = sqrt(r2);

    return (r - fh->R1)*fh->RS + fh->Smin;                      // Linear interpolation from Smin to Smax
//  return tanh(((r - hhmd->R1)*hhmd->RS - 0.5)*6.0)*0.5 + 0.5; // Smoother interpolation
}


double fhmd_Sxyz_d(const dvec x, const dvec c, FHMD *fh)
{
    rvec xr;
    copy_dvec_to_rvec(x, xr);
    return fhmd_Sxyz_r(xr, c, fh);
}


/*
 ******************** Estimate S in the cells and cell faces ********************
 */
void FH_S_weighted(FHMD *fh)
{
    for(int i = 0; i < fh->Ntot; i++)
    {
        if(fh->S_function != constant_S)
            fh->arr[i].S = 1 - fh->arr[i].ro_md_s/fh->arr[i].ro_md;
        else
            fh->arr[i].S = fh->S;
    }

    ivec ind;

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                    fh->arr[L].Sf[d] = (fh->arr[CL].S + fh->arr[C].S)*0.5;
            }
        }
    }
}


void FH_S_precise(FHMD *fh)
{
    dvec coordl, coordr;
    ivec ind;

    for(int k = fh->N_shift[2]; k < (fh->N_md[2] + fh->N_shift[2]); k++)
    {
        for(int j = fh->N_shift[1]; j < (fh->N_md[1] + fh->N_shift[1]); j++)
        {
            for(int i = fh->N_shift[0]; i < (fh->N_md[0] + fh->N_shift[0]); i++)
            {
                ASSIGN_IND(ind, i, j, k);

                if(fh->S_function == moving_sphere)
                    fh->arr[C].S = fhmd_Sxyz_d(fh->grid.c[C], fh->protein_com, fh);         // MD/FH sphere follows protein
                else if(fh->S_function == fixed_sphere)
                    fh->arr[C].S = fhmd_Sxyz_d(fh->grid.c[C], fh->box05, fh);               // Fixed MD/FH sphere
                else
                    fh->arr[C].S = fh->S;                                                   // Constant S

                for(int d = 0; d < DIM; d++)
                {
                    if(fh->S_function != constant_S)
                    {
                        switch(d)
                        {
                        case 0:
                            ASSIGN_DVEC(coordl, fh->grid.n[C][0], fh->grid.c[C][1], fh->grid.c[C][2]);
                            ASSIGN_DVEC(coordr, fh->grid.n[CR][0], fh->grid.c[C][1], fh->grid.c[C][2]);
                            break;
                        case 1:
                            ASSIGN_DVEC(coordl, fh->grid.c[C][0], fh->grid.n[C][1], fh->grid.c[C][2]);
                            ASSIGN_DVEC(coordr, fh->grid.c[C][0], fh->grid.n[CR][1], fh->grid.c[C][2]);
                            break;
                        case 2:
                            ASSIGN_DVEC(coordl, fh->grid.c[C][0], fh->grid.c[C][1], fh->grid.n[C][2]);
                            ASSIGN_DVEC(coordr, fh->grid.c[C][0], fh->grid.c[C][1], fh->grid.n[CR][2]);
                            break;
                        }

                        if(fh->S_function == moving_sphere)
                        {
                            fh->arr[L].Sf[d] = fhmd_Sxyz_d(coordl, fh->protein_com, fh);    // MD/FH sphere follows protein
                            fh->arr[R].Sf[d] = fhmd_Sxyz_d(coordr, fh->protein_com, fh);
                        }
                        else
                        {
                            fh->arr[L].Sf[d] = fhmd_Sxyz_d(coordl, fh->box05, fh);          // Fixed MD/FH sphere
                            fh->arr[R].Sf[d] = fhmd_Sxyz_d(coordr, fh->box05, fh);
                        }
                    }
                    else
                    {
                            fh->arr[L].Sf[d] = fh->S;                                       // Constant S
                            fh->arr[R].Sf[d] = fh->S;
                    }
                }

            }
        }
    }
}


void FH_S(FHMD *fh)
{
    for(int i = 0; i < fh->Ntot; i++)
    {
        if(fh->S_function == moving_sphere)
            fh->arr[i].S = fhmd_Sxyz_d(fh->grid.c[i], fh->protein_com, fh);         // MD/FH sphere follows protein
        else if(fh->S_function == fixed_sphere)
            fh->arr[i].S = fhmd_Sxyz_d(fh->grid.c[i], fh->box05, fh);               // Fixed MD/FH sphere
        else
            fh->arr[i].S = fh->S;                                                   // Constant S

        if(fh->grid.md[i] == FH_zone)
            fh->arr[i].S = 1;                                                       // Pure FH zone
    }

    ivec ind;

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                    fh->arr[L].Sf[d] = (fh->arr[CL].S + fh->arr[C].S)*0.5;
            }
        }
    }
/*
    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                if(fh->grid.md[C] == FH_zone)
                {
                    for(int d = 0; d < DIM; d++)
                    {
                        fh->arr[L].Sf[d] = 1;
                        fh->arr[R].Sf[d] = 1;
                    }
                }
            }
        }
    }
*/
}


/*
 ******************** Find protein molecule in the box ********************
 */
void fhmd_find_protein(gmx_mtop_t *mtop, int N_atoms, real mass[], t_commrec *cr, FHMD *fh)
{
    int    ind, res_nr;
    char  *atomname, *resname;
    int    protein_n    = 0;
    double protein_mass = 0;

    for(int n = 0; n < N_atoms; n++)
    {
        if(PAR(cr) && DOMAINDECOMP(cr))
        {
            ind = cr->dd->gatindex[n];
        } else {
            ind = n;
        }

        gmx_mtop_atominfo_global(mtop, ind, &atomname, &res_nr, &resname);

        if(strcmp(resname, "SOL") && strcmp(resname, "CL") && strcmp(resname, "NA"))
        {
            protein_n++;
            protein_mass += mass[n];
        }
    }

    if(PAR(cr))
    {
        gmx_sumi(1, &protein_n, cr);
        gmx_sumd(1, &protein_mass, cr);
    }

    fh->protein_N    = protein_n;
    fh->protein_mass = protein_mass;
}


/*
 ******************** Find protein molecule centre of mass ********************
 */
void fhmd_find_protein_com(gmx_mtop_t *mtop, int N_atoms, rvec x[], real mass[], t_commrec *cr, FHMD *fh)
{
    int    ind, res_nr;
    char  *atomname, *resname;
    rvec   pcom;
    dvec   r, rm;
    double xd;

    clear_dvec(rm);

    for(int n = 0; n < N_atoms; n++)
    {
        if(PAR(cr) && DOMAINDECOMP(cr))
        {
            ind = cr->dd->gatindex[n];
        } else {
            ind = n;
        }

        gmx_mtop_atominfo_global(mtop, ind, &atomname, &res_nr, &resname);

        if(strcmp(resname, "SOL") && strcmp(resname, "CL") && strcmp(resname, "NA"))
        {
            for(int d = 0; d < DIM; d++)
            {
                pcom[d] = fh->protein_com[d];

                xd = x[n][d] - pcom[d];

                if(fabs(xd) <= fh->box05[d])
                    r[d] = x[n][d];
                else if(xd > fh->box05[d])
                    r[d] = x[n][d] - fh->box[d];
                else    // In case if(xd < -fh->box05[d])
                    r[d] = x[n][d] + fh->box[d];

                rm[d] += mass[n]*r[d];
            }
        }
    }

    if(PAR(cr))
    {
        gmx_sumd(3, rm, cr);
    }

    for(int d = 0; d < DIM; d++)
        pcom[d] = rm[d]/fh->protein_mass;

    PBC(fh->protein_com, pcom, fh->box);

#ifdef FHMD_DEBUG_COM
    if(MASTER(cr) && !(fh->step_MD % 10000))
        printf("FHMD DEBUG: Protein COM position: %g, %g, %g\n", fh->protein_com[0], fh->protein_com[1], fh->protein_com[2]);
#endif
}
