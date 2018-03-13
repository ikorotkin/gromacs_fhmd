#ifndef FHMD_ESTIMATE_H_
#define FHMD_ESTIMATE_H_

void fhmd_reset_statistics(FHMD *fh);
void fhmd_collect_statistics(FHMD *fh);
void fhmd_update_statistics(FHMD *fh);
void fhmd_print_statistics(FHMD *fh, t_commrec *cr);

void fhmd_couette_avg(FHMD *fh, rvec x[], rvec v[], int N_atoms, t_commrec *cr);

#endif /* FHMD_ESTIMATE_H_ */
