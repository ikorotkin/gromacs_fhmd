#ifndef FHMD_COUPLING_H_
#define FHMD_COUPLING_H_

void fhmd_update_MD_in_FH(rvec x[], rvec v[], real mass[], rvec f[], int N_atoms, FHMD *fh);
void fhmd_sum_arrays(t_commrec *cr, FHMD *fh);
void fhmd_calculate_MDFH_terms(FHMD *fh);

#endif /* FHMD_COUPLING_H_ */
