#ifndef FHMD_INIT_H_
#define FHMD_INIT_H_

int fhmd_init(matrix box, int N_atoms, real mass[], rvec x[], double dt_md, gmx_mtop_t *mtop, t_commrec *cr, FHMD *fh);

#endif /* FHMD_INIT_H_ */
