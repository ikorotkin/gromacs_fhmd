#ifndef FHMD_SFUNCTION_H_
#define FHMD_SFUNCTION_H_

double fhmd_Sxyz_r(const rvec x, const dvec c, FHMD *fh);
double fhmd_Sxyz_d(const dvec x, const dvec c, FHMD *fh);

void FH_S(FHMD *fh);
void FH_S_weighted(FHMD *fh);
void FH_S_precise(FHMD *fh);

void fhmd_find_protein(gmx_mtop_t *mtop, int N_atoms, real mass[], t_commrec *cr, FHMD *fh);
void fhmd_find_protein_com(gmx_mtop_t *mtop, int N_atoms, rvec x[], real mass[], t_commrec *cr, FHMD *fh);

#endif /* FHMD_SFUNCTION_H_ */
