#ifndef FHMD_SFUNCTION_H_
#define FHMD_SFUNCTION_H_

double fhmd_Sxyz_r(const rvec x, const dvec c, FHMD *fh);

void fhmd_find_protein_com(gmx_mtop_t *mtop, int N_atoms, rvec x[], real mass[], t_commrec *cr, FHMD *fh);

#endif /* FHMD_SFUNCTION_H_ */
