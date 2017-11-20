#ifndef FHMD_PARSER_H_
#define FHMD_PARSER_H_

int  parse_prm(char const *fname, FHMD *fh);
void skip_line(FILE *fprm);
int  assign_int_value(int *v, char *line, FILE *fprm);
int  assign_double_value(double *v, char *line, FILE *fprm);

#endif /* FHMD_PARSER_H_ */
