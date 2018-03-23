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
/*    write_dump(f_fh[0],       "f_fh_X");
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

    write_dump(ro_prime,      "ro_prime");
    write_dump(ro_star,       "ro_star");
    write_dump(m_prime[0],    "m_prime_X");
    write_dump(m_prime[1],    "m_prime_Y");
    write_dump(m_prime[2],    "m_prime_Z");
    write_dump(m_star[0],     "m_star_X");
    write_dump(m_star[1],     "m_star_Y");
    write_dump(m_star[2],     "m_star_Z");

    write_dump(S,             "S");
    write_dump(Sf[0],         "Sf_X");
    write_dump(Sf[1],         "Sf_Y");
    write_dump(Sf[2],         "Sf_Z");

    write_dump(ux_avg,        "Ux_tilde_avg");
    write_dump(mx_avg,        "mx_star_avg");
*/
#endif
}


void fhmd_couette_avg_write(FHMD *fh)
{
    FILE *fw;

    const double eps = 1e-15;

    if(fh->step_MD == 0)
    {
        fw = fopen("couette.txt", "w");
    } else {
        fw = fopen("couette.txt", "a");
    }

    fprintf(fw, "MD step = %d\n", fh->step_MD);
    fprintf(fw, "layer\ty, nm\t<Ux>, nm/ps\t<Ux> (MD), nm/ps\tUx, nm/ps\tUx (MD), nm/ps\n");

    for(int k = 0; k < FHMD_COUETTE_LAYERS; k++)
    {
        fprintf(fw, "%d\t%g\t%g\t%g\t%g\t%g\n", k, ((double)(k) + 0.5)/(double)(FHMD_COUETTE_LAYERS)*fh->box[YY],
                fh->avg_vel_tot[k]/fh->avg_n_tot[k], fh->avg_vel_S_tot[k]/fh->avg_n_tot[k],
                fh->avg_vel[k]/((double)(fh->avg_n[k]) + eps), fh->avg_vel_S[k]/((double)(fh->avg_n_S[k]) + eps));
    }

    fprintf(fw, "\n");
    fclose(fw);
}


void fhmd_acoustic_avg_write(FHMD *fh)
{
    FILE *fw;
    double invn;

    if(fh->step_MD == 0)
    {
        fw = fopen("acoustic.txt", "w");
    } else {
        fw = fopen("acoustic.txt", "a");
    }

    fprintf(fw, "MD step = %d\n", fh->step_MD);
    fprintf(fw, "cell\tx, nm\t<U_tilde>, nm/ps\t<Rho_tilde>, nm/ps\t<m_star/rho>, nm/ps\n");

    for(int i = 0; i < NX; i++)
    {
        fprintf(fw, "%d\t%g\t%g\t%g\t%g\n", i, ((double)(i) + 0.5)*fh->Lx/(double)(NX),
                fh->arr[i].ux_avg, fh->arr[i].ro_avg, fh->arr[i].mx_avg);
    }

    fprintf(fw, "\n");
    fclose(fw);

    if(fh->step_MD == 0)
    {
        fw = fopen("acoustic_t.txt", "w");
    } else {
        fw = fopen("acoustic_t.txt", "a");
    }

    fprintf(fw, "MD step = %d\n", fh->step_MD);
    fprintf(fw, "t, ps\tUx1\tUx2\tUx3\tUx4\tUx5\tRo1\tRo2\tRo3\tRo4\tRo5\n");

    for(int j = 0; j < FHMD_COUETTE_LAYERS; j++)
    {
        if(fh->ac_n[j] > 0)
            invn = 1.0/(double)(fh->ac_n[j]);
        else
            invn = 0;

        fprintf(fw, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", (double)(j)/(double)(FHMD_COUETTE_LAYERS)*556.61,
                fh->ac_vel[0][j]*invn, fh->ac_vel[1][j]*invn, fh->ac_vel[2][j]*invn, fh->ac_vel[3][j]*invn, fh->ac_vel[4][j]*invn,
                fh->ac_rho[0][j]*invn, fh->ac_rho[1][j]*invn, fh->ac_rho[2][j]*invn, fh->ac_rho[3][j]*invn, fh->ac_rho[4][j]*invn);
    }

    fprintf(fw, "\n");
    fclose(fw);
}


void writePoint(FILE *fout, FHMD *fh, const char *line, char ch)
{
    int c = 0;
    ivec ind;

    fprintf(fout, "%s", line);

    for(int k = 0; k <= NZ; k++) {
        for(int j = 0; j <= NY; j++) {
            for(int i = 0; i <= NX; i++) {
                ASSIGN_IND(ind, i - (int)(i/NX), j - (int)(j/NY), k - (int)(k/NZ));
                     if(ch == 'X' && i < NX)  fprintf(fout, "%e ", fh->grid.n[C][0]);
                else if(ch == 'X' && i == NX) fprintf(fout, "%e ", fh->grid.n[C][0] + fh->grid.h[C][0]);
                else if(ch == 'Y' && j < NY)  fprintf(fout, "%e ", fh->grid.n[C][1]);
                else if(ch == 'Y' && j == NY) fprintf(fout, "%e ", fh->grid.n[C][1] + fh->grid.h[C][1]);
                else if(ch == 'Z' && k < NZ)  fprintf(fout, "%e ", fh->grid.n[C][2]);
                else if(ch == 'Z' && k == NZ) fprintf(fout, "%e ", fh->grid.n[C][2] + fh->grid.h[C][2]);
                if(c++ == 10) {fprintf(fout, "\n"); c = 0;}
            }
        }
    }
}


void writeCons(FILE *fout, FHMD *fh, const char *line, char ch)
{
    int c = 0;
    ivec ind;

    fprintf(fout, "%s", line);

    for(int k = 0; k < NZ; k++) {
        for(int j = 0; j < NY; j++) {
            for(int i = 0; i < NX; i++) {
                ASSIGN_IND(ind, i, j, k);
                     if(ch == 'U') fprintf(fout, "%e ", fh->arr[C].ux_avg); // .u_fh[0]); !!!
                else if(ch == 'V') fprintf(fout, "%e ", fh->arr[C].u_fh[1]);
                else if(ch == 'W') fprintf(fout, "%e ", fh->arr[C].u_fh[2]);
                else if(ch == 'R') fprintf(fout, "%e ", fh->arr[C].ro_fh);
                if(c++ == 10) {fprintf(fout, "\n"); c = 0;}
            }
        }
    }
}


// TODO: This function is a copy/paste from the previous code -- should be rewritten completely
void fhmd_write_tecplot_data(FHMD *fh, int step, double time)
{
    FILE *fout;                                                     // Tecplot data file

    static int zc = 0;

    char fname[32];

    sprintf(fname, "tecplot/data_%7.7d.dat", step);                 // Tecplot data filename

    if((fout = fopen(fname, "w")) == NULL)                          // Open file
        printf("\n ERROR creating %s for output!\n", fname);

    zc++;

    fprintf(fout, "TITLE=\"OUT\"\nVARIABLES=\"X\",\"Y\",\"Z\",\"U_tilde\",\"V_tilde\",\"W_tilde\",\"RHO_tilde\"\n");
    fprintf(fout, "ZONE T=\"%f\", I=%d, J=%d, K=%d, DATAPACKING=BLOCK\nVARLOCATION=([4-7]=CELLCENTERED)\n", time, NX+1, NY+1, NZ+1);
    fprintf(fout, "STRANDID=%d\nSOLUTIONTIME=%g", zc, time);

    writePoint(fout, fh, "\n\n# X:\n", 'X');                        // Write X
    writePoint(fout, fh, "\n\n# Y:\n", 'Y');                        // Write Y
    writePoint(fout, fh, "\n\n# Z:\n", 'Z');                        // Write Z

    writeCons(fout, fh, "\n\n# U:\n",   'U');                       // Write U
    writeCons(fout, fh, "\n\n# V:\n",   'V');                       // Write V
    writeCons(fout, fh, "\n\n# W:\n",   'W');                       // Write W
    writeCons(fout, fh, "\n\n# RHO:\n", 'R');                       // Write RHO

    fclose(fout);
}


#ifdef FHMD_PARAVIEW
void write_paraview_data(HHMD *hhmd) {

    const int ln = 5;       // Numbers per one line

    FILE *fout;
    char fname[32];
    static int zc = 0;
    int lc = 1;

    int nx = hhmd->cellNum.x;
    int ny = hhmd->cellNum.y;
    int nz = hhmd->cellNum.z;

    sprintf(fname, "paraview_%d.vtk", zc);

    if((fout = fopen(fname, "w")) == NULL)
        printf("\n ERROR creating %s for output!\n", fname);

    zc++;

    fprintf(fout, "# vtk DataFile Version 3.0\nData Output\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS %d %d %d\n", nx + 1, ny + 1, nz + 1);
    fprintf(fout, "POINTS %d float\n", (nx + 1)*(ny + 1)*(nz + 1));

    for(int k = 0; k < nz + 1; k++) {
        for(int j = 0; j < ny + 1; j++) {
            for(int i = 0; i < nx + 1; i++) {
                fprintf(fout, "%f %f %f\n", hhmd->grid.n[i].x, hhmd->grid.n[j].y, hhmd->grid.n[k].z);
            }
        }
    }

    fprintf(fout, "CELL_DATA %d\n", nx*ny*nz);
    fprintf(fout, "FIELD FieldData %d\n", 1);
    fprintf(fout, "Density %d %d float\n", 1, nx*ny*nz);

    for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
            for(int i = 0; i < nx; i++) {
                fprintf(fout, "%f", hhmd->arr->ro_fh[i][j][k]);
                if(lc++ >= ln) {
                    fprintf(fout, "\n");
                    lc = 1;
                } else {
                    fprintf(fout, " ");
                }
            }
        }
    }

    fprintf(fout, "VECTORS Velocity float\n");

    lc = 1;

    for(int k = 0; k < nz; k++) {
        for(int j = 0; j < ny; j++) {
            for(int i = 0; i < nx; i++) {
                fprintf(fout, "%f %f %f", hhmd->arr->u_fh[i][j][k].x, hhmd->arr->u_fh[i][j][k].y, hhmd->arr->u_fh[i][j][k].z);
                if(lc++ >= ln) {
                    fprintf(fout, "\n");
                    lc = 1;
                } else {
                    fprintf(fout, "   ");
                }
            }
        }
    }

    fprintf(fout, "POINT_DATA %d\n", (nx + 1)*(ny + 1)*(nz + 1));

    fclose(fout);

}
#endif

