#include "data_structures.h"
#include "macro.h"

#include "gromacs/topology/mtop_util.h"     /* This is for gmx_mtop_atominfo_global() */


double fhmd_Sxyz_r(const rvec x, const dvec c, FHMD *fh)
{
    rvec xd;

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

    ZERO(rm);

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
    if(MASTER(cr) && !(fh->step_MD % 100))
        printf("FHMD DEBUG: Protein COM position: %g, %g, %g\n", fh->protein_com[0], fh->protein_com[1], fh->protein_com[2]);
#endif
}
