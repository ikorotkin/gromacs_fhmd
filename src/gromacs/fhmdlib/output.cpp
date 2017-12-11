#include "data_structures.h"
#include "macro.h"


#define write_dump(var, name) \
    sprintf(fname, "dump_%s.txt", name); \
    if(fh->step_MD == 0) { \
        fw = fopen(fname, "w"); \
        write_header(fw, fh); \
    } else { \
        fw = fopen(fname, "a"); \
    } \
    fprintf(fw, "\n%d\t", fh->step_MD); \
    for(int i = 0; i < fh->Ntot; i++) fprintf(fw, "%g\t", fh->arr[i].var); \
    fclose(fw);


void write_header(FILE *fw, FHMD *fh)
{
    fprintf(fw, "step\t");

    for(int k = 0; k < NZ; k++)
        for(int j = 0; j < NY; j++)
            for(int i = 0; i < NX; i++)
                fprintf(fw, "cell %d-%d-%d\t", i, j, k);
}


void fhmd_dump_all(FHMD *fh)
{
    FILE *fw;
    char  fname[64];

    write_dump(ro_md,   "ro_md");
    write_dump(ro_fh,   "ro_fh");
    write_dump(u_md[0], "u_md_X");
    write_dump(u_md[1], "u_md_Y");
    write_dump(u_md[2], "u_md_Z");
    write_dump(u_fh[0], "u_fh_X");
    write_dump(u_fh[1], "u_fh_Y");
    write_dump(u_fh[2], "u_fh_Z");

#ifdef FHMD_DEBUG
    write_dump(f_fh[0],       "f_fh_X");
    write_dump(f_fh[1],       "f_fh_Y");
    write_dump(f_fh[2],       "f_fh_Z");
    write_dump(alpha_term[0], "alpha_term_X");
    write_dump(alpha_term[1], "alpha_term_Y");
    write_dump(alpha_term[2], "alpha_term_Z");
    write_dump(beta_term[0],  "beta_term_X");
    write_dump(beta_term[1],  "beta_term_Y");
    write_dump(beta_term[2],  "beta_term_Z");
    write_dump(grad_ro[0],    "grad_ro_X");
    write_dump(grad_ro[1],    "grad_ro_Y");
    write_dump(grad_ro[2],    "grad_ro_Z");
#endif
}
