#ifndef FHMD_FH_FUNCTIONS_H_
#define FHMD_FH_FUNCTIONS_H_

void FH_init(FHMD *fh, t_commrec *cr);
void FH_predictor(FHMD *fh);
void FH_corrector(FHMD *fh);
void FH_char(FHMD *fh);
void FH_do_single_timestep(FHMD *fh);
void FH_equilibrate(FHMD *fh);
void define_FH_grid(t_commrec *cr, FHMD *fh);

#endif /* FHMD_FH_FUNCTIONS_H_ */
