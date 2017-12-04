#include "data_structures.h"
#include "parser.h"


int parse_prm(char const *fname, FHMD *fh)
{
    const char c = ';';     // comment delimiter

    FILE *fprm;

    if((fprm = fopen(fname, "r")) == NULL) return -1;   // file not found

    char line[255];
    int  ok = 1;

    while(fscanf(fprm, "%s", line) != -1)
    {
        if(line[0] == c) skip_line(fprm);   // skip comment

        if(!strcmp(line, "S") || !strcmp(line, "s"))
            ok = assign_double_value(&fh->S, line, fprm);
        else if(!strcmp(line, "R1"))
            ok = assign_double_value(&fh->R1, line, fprm);
        else if(!strcmp(line, "R2"))
            ok = assign_double_value(&fh->R2, line, fprm);
        else if(!strcmp(line, "Smin"))
            ok = assign_double_value(&fh->Smin, line, fprm);
        else if(!strcmp(line, "Smax"))
            ok = assign_double_value(&fh->Smax, line, fprm);
        else if(!strcmp(line, "alpha"))
            ok = assign_double_value(&fh->alpha, line, fprm);
        else if(!strcmp(line, "beta"))
            ok = assign_double_value(&fh->beta, line, fprm);
        else if(!strcmp(line, "Nx"))
            ok = assign_int_value(&fh->N[0], line, fprm);
        else if(!strcmp(line, "Ny"))
            ok = assign_int_value(&fh->N[1], line, fprm);
        else if(!strcmp(line, "Nz"))
            ok = assign_int_value(&fh->N[2], line, fprm);
        else if(!strcmp(line, "FH_EOS"))
            ok = assign_int_value(&fh->FH_EOS, line, fprm);
        else if(!strcmp(line, "FH_equil"))
            ok = assign_int_value(&fh->FH_equil, line, fprm);
        else if(!strcmp(line, "FH_step"))
            ok = assign_int_value(&fh->FH_step, line, fprm);
        else if(!strcmp(line, "FH_dens"))
            ok = assign_double_value(&fh->FH_dens, line, fprm);
        else if(!strcmp(line, "FH_temp"))
            ok = assign_double_value(&fh->FH_temp, line, fprm);

        if(!ok) return 0;   // error in prm-file
    }

    fclose(fprm);

    return 1;
}


void skip_line(FILE *fprm)
{
    int ok = 1;
    char c = '0';

    while((c != '\n') && (ok == 1)) ok = fscanf(fprm, "%c", &c);
}


int assign_int_value(int *v, char *line, FILE *fprm)
{
    int ok = fscanf(fprm, "%s", line);  // skip 'equal' delimiter

    ok = fscanf(fprm, "%s", line);

    if(sscanf(line, "%d", v) != 1) {
        return 0;
    } else {
        return 1;
    }
}


int assign_double_value(double *v, char *line, FILE *fprm)
{
    float f;

    int ok = fscanf(fprm, "%s", line);  // skip 'equal' delimiter

    ok = fscanf(fprm, "%s", line);

    if(sscanf(line, "%f", &f) != 1) {
        return 0;
    } else {
        *v = f;
        return 1;
    }
}
