#ifndef FHMD_NEW_UPDATE_H_
#define FHMD_NEW_UPDATE_H_

void fhmd_update_coords(FILE        *fplog,
                   gmx_int64_t       step,
                   t_inputrec       *inputrec,  /* input record and box stuff   */
                   t_mdatoms        *md,
                   t_state          *state,
                   rvec             *f,    /* forces on home particles */
                   t_fcdata         *fcd,
                   gmx_ekindata_t   *ekind,
                   matrix            M,
                   gmx_update_t     *upd,
                   int               UpdatePart,
                   t_commrec        *cr, /* these shouldn't be here -- need to think about it */
                   gmx_constr_t      constr,
                   FHMD             *fh);

void fhmd_do_update_md(int start, int nrend,
                       double dt, int nstpcouple,
                       t_grp_tcstat *tcstat,
                       double nh_vxi[],
                       gmx_bool bNEMD, t_grp_acc *gstat, rvec accel[],
                       ivec nFreeze[],
                       real invmass[],
                       unsigned short ptype[], unsigned short cFREEZE[],
                       unsigned short cACC[], unsigned short cTC[],
                       rvec x[], rvec xprime[], rvec v[],
                       rvec f[], matrix M,
                       gmx_bool bNH, gmx_bool bPR, t_commrec *cr, FHMD *fh);

#endif /* FHMD_NEW_UPDATE_H_ */
