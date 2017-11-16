/*

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int assign_int_value(int *v, char *line, FILE *fprm) {

    fscanf(fprm, "%s", line);   // skip equal delimeter
    fscanf(fprm, "%s", line);
    if(sscanf(line, "%d", v) != 1) return 0;
    else return 1;

}


int assign_double_value(double *v, char *line, FILE *fprm) {
    
    float f;
    
    fscanf(fprm, "%s", line);   // skip equal delimeter
    fscanf(fprm, "%s", line);
    if(sscanf(line, "%f", &f) != 1) {
        return 0;
    } else {
        *v = f;
        return 1;
    }

}


int parser(char *fname, hmd_global *prm) {
    
    const char c = ';';     // comment delimeter

    FILE *fprm;

    if((fprm = fopen("coupling.txt", "r")) == NULL) return 1;

    char line[255];
    int read_value = 0, v = 0, argon = 0, ok;
    double d = 0;

    while(fscanf(fprm, "%s", line) != -1)
    {

        if(line[0] == c) continue;  // skip comment
        if(line[0] == e) {read_value = 1; continue;}

        if(!strcmp(line, "argon")) ok=assign_int_value(&argon, line, fprm);
        else if(!strcmp(line, "water")) ok=assign_double_value(&d, line, fprm);

        if(!ok) printf("ERROR");

        printf("line: %s\n", line);

            printf("v=%d\n", argon);
            printf("w=%g\n", d);

    }
    
    return 1;
    
}


int main(int argc, char **argv) {

    int ok = parser("coupling.prm", prm);

}
*/

