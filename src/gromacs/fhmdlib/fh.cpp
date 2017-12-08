#include "data_structures.h"
#include "macro.h"
#include "fh.h"


void FH_init(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    switch(fh->eos)
    {
    case eos_argon:
        MU     = 0;
        KAPPA  = 0;
        EOS_A  = 0;
        EOS_B  = 0;
        EOS_C  = 0;
        EOS_D  = 0;
        EOS_E  = 0;
        P_INIT = 0;
        SOUND  = 0;
        break;
    case eos_spce:
        MU     = FHMD_EOS_SPCE_MU;
        KAPPA  = FHMD_EOS_SPCE_KAPPA;
        EOS_A  = FHMD_EOS_SPCE_A;
        EOS_B  = FHMD_EOS_SPCE_B;
        EOS_C  = FHMD_EOS_SPCE_C;
        P_INIT = fh->FH_dens*(fh->FH_dens*EOS_A + EOS_B) + EOS_C;
        SOUND  = sqrt(2.0*fh->FH_dens*EOS_A + EOS_B);
        break;
    }

#ifdef FHMD_DEBUG
    printf(MAKE_YELLOW "FHMD DEBUG: MU = %g, KAPPA = %g\n", MU, KAPPA);
    printf(MAKE_YELLOW "FHMD DEBUG: A = %g, B = %g, C = %g\n", EOS_A, EOS_B, EOS_C);
    printf(MAKE_YELLOW "FHMD DEBUG: Initial Pressure = %g, Speed of Sound = %g\n", P_INIT, SOUND*1000.0);
    printf(RESET_COLOR "\n");
#endif

    // Viscosity terms for viscous stress
    VISC1 = 4.0/3.0 + KAPPA/MU;
    VISC2 = KAPPA/MU - 2.0/3.0;

    for(int i = 0; i < fh->Ntot; i++)
    {
        arr[i].ro_fh   = fh->FH_dens;
        arr[i].ro_fh_n = fh->FH_dens;
        arr[i].p       = P_INIT;
        arr[i].pn      = P_INIT;

        for(int d = 0; d < DIM; d++)
        {
            arr[i].rof[d]  = fh->FH_dens;
            arr[i].rofn[d] = fh->FH_dens;
            arr[i].pf[d]   = P_INIT;
            arr[i].pfn[d]  = P_INIT;
        }
    }

    T        = 0;       // Current time
    STEP     = 0;       // Current time step

    if(fh->FH_blend >= 0)
        blend = fh->FH_blend;

    // Reset statistics
    T_AVG    = 0;
    NT_AVG   = 0;
    RHO_AVG  = 0;
    RHO2_AVG = 0;
    P_AVG    = 0;
    N_AVG    = 0;

    for(int d = 0; d < DIM; d++)
    {
        U_AVG[d]  = 0;
        U2_AVG[d] = 0;
    }
}


void FH_predictor(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    ivec   ind;
    double DT = fh->dt_FH;
    dvec   F, PG, TAU, TAURAN;
    matrix TAUL, TAUR;

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                {
                    // Mass flux
                    F[d] = (arr[R].rof[d]*arr[R].uf[d][d] - arr[L].rof[d]*arr[L].uf[d][d])/fh->grid.h[C][d];
                }

                // Mass conservation
                arr[C].ro_fh_n = arr[C].ro_fh - 0.5*DT*SUM(F);

                for(int d = 0; d < DIM; d++)
                {
                    // Pressure gradient
                    PG[d] = (arr[R].pf[d] - arr[L].pf[d])/fh->grid.h[C][d];
                }

                for(int dim = 0; dim < DIM; dim++)
                {
                    VISCOUS_FLUX(u_fh);     // UNIFORM GRID!

                    // Momentum flux, viscous and random stress
                    for(int d = 0; d < DIM; d++)
                    {
                        F[d]      = (arr[R].rof[d]*arr[R].uf[d][d]*arr[R].uf[dim][d] - arr[L].rof[d]*arr[L].uf[d][d]*arr[L].uf[dim][d])/fh->grid.h[C][d];
                        TAU[d]    = (TAUR[dim][d] - TAUL[dim][d])/fh->grid.h[C][d];
                        TAURAN[d] = (arr[R].rans[dim][d] - arr[L].rans[dim][d])/fh->grid.h[C][d];
                    }

                    // FH force
                    arr[C].f_fh[dim] = -SUM(F) - PG[dim] + SUM(TAU) + SUM(TAURAN);

                    // Momentum conservation
                    arr[C].u_fh_n[dim] = (arr[C].ro_fh*arr[C].u_fh[dim] + 0.5*DT*arr[C].f_fh[dim])/arr[C].ro_fh_n;
                }

                // Pressure
                switch(fh->eos)
                {
                case eos_argon:
                    arr[C].pn = EOS_D*(1./3.*EOS_A*EOS_A*pow(EOS_E*arr[C].ro_fh_n - EOS_B, 3) + EOS_C);
                    break;
                case eos_spce:
                    arr[C].pn = arr[C].ro_fh_n*(arr[C].ro_fh_n*EOS_A + EOS_B) + EOS_C;
                    break;
                }
            }
        }
    }
}


void FH_corrector(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    ivec   ind;
    double DT = fh->dt_FH;
    dvec   F, PG, TAU, TAURAN;
    matrix TAUL, TAUR;

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                {
                    // Mass flux
                    F[d] = (arr[R].rofn[d]*arr[R].ufn[d][d] - arr[L].rofn[d]*arr[L].ufn[d][d])/fh->grid.h[C][d];
                }

                // Mass conservation
                arr[C].ro_fh = arr[C].ro_fh_n - 0.5*DT*SUM(F);

                for(int d = 0; d < DIM; d++)
                {
                    // Pressure gradient
                    PG[d] = (arr[R].pfn[d] - arr[L].pfn[d])/fh->grid.h[C][d];
                }

                for(int dim = 0; dim < DIM; dim++)
                {
                    VISCOUS_FLUX(u_fh_n);     // UNIFORM GRID!

                    // Momentum flux, viscous and random stress
                    for(int d = 0; d < DIM; d++)
                    {
                        F[d]      = (arr[R].rofn[d]*arr[R].ufn[d][d]*arr[R].ufn[dim][d] - arr[L].rofn[d]*arr[L].ufn[d][d]*arr[L].ufn[dim][d])/fh->grid.h[C][d];
                        TAU[d]    = (TAUR[dim][d] - TAUL[dim][d])/fh->grid.h[C][d];
                        TAURAN[d] = (arr[R].rans[dim][d] - arr[L].rans[dim][d])/fh->grid.h[C][d];
                    }

                    // FH force
                    arr[C].f_fh[dim] = -SUM(F) - PG[dim] + SUM(TAU) + SUM(TAURAN);

                    // Momentum conservation
                    arr[C].u_fh[dim] = (arr[C].ro_fh_n*arr[C].u_fh_n[dim] + 0.5*DT*arr[C].f_fh[dim])/arr[C].ro_fh;
                }

                // Pressure
                switch(fh->eos)
                {
                case eos_argon:
                    arr[C].p = EOS_D*(1./3.*EOS_A*EOS_A*pow(EOS_E*arr[C].ro_fh - EOS_B, 3) + EOS_C);
                    break;
                case eos_spce:
                    arr[C].p = arr[C].ro_fh*(arr[C].ro_fh*EOS_A + EOS_B) + EOS_C;
                    break;
                }
            }
        }
    }
}


void FH_char(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    const double RO1 = 1.0/fh->FH_dens;
    const double DT  = fh->dt_FH;

    ivec ind, indR, indL;
    int  d1, d2;

    double RCLN, RCL, QCLN, QCL, VCLN, VCL, VL, VPL, WCLN, WCL, WL, WPL, RL, RPL, QL, QPL;
    double P1, P2, P3, RCRN, QCRN, RCR, QCR, RR, QR, VCRN, VCR, VR, WCRN, WCR, WR, RPR, QPR, VPR, WPR;
    double L1L, L2L, L3L, L1R, L2R, L3R, RPNL, QPNL, VPNL, WPNL, RPNR, QPNR, VPNR, WPNR;
    double QQ1, QQ2, QQ3, QQ4, RMAX, QMAX, VMAX, WMAX, RMIN, QMIN, VMIN, WMIN;
    double RPN1, QPN1, VPN1, WPN1, RPN2, QPN2, VPN2, WPN2, RPN, QPN, VPN, WPN;

    for(int d = 0; d < DIM; d++)
    {
        d1 = d  + 1; if(d1 >= DIM) d1 = 0;
        d2 = d1 + 1; if(d2 >= DIM) d2 = 0;

        for(int k = 0; k < fh->N[d2]; k++)
        {
            for(int j = 0; j < fh->N[d1]; j++)
            {

                RCLN = arr[L0].u_fh_n[d] + SOUND*log(arr[L0].ro_fh_n*RO1);
                RCL  = arr[L0].u_fh[d]   + SOUND*log(arr[L0].ro_fh*RO1);

                QCLN = arr[L0].u_fh_n[d] - SOUND*log(arr[L0].ro_fh_n*RO1);
                QCL  = arr[L0].u_fh[d]   - SOUND*log(arr[L0].ro_fh*RO1);

                VCLN = arr[L0].u_fh_n[d1];
                VCL  = arr[L0].u_fh[d1];
                VL   = arr[L0].uf[d1][d];
                VPL  = arr[L1].uf[d1][d];

                WCLN = arr[L0].u_fh_n[d2];
                WCL  = arr[L0].u_fh[d2];
                WL   = arr[L0].uf[d2][d];
                WPL  = arr[L1].uf[d2][d];

                RL   = arr[L0].uf[d][d] + SOUND*log(arr[L0].rof[d]*RO1);
                RPL  = arr[L1].uf[d][d] + SOUND*log(arr[L1].rof[d]*RO1);
                QL   = arr[L0].uf[d][d] - SOUND*log(arr[L0].rof[d]*RO1);
                QPL  = arr[L1].uf[d][d] - SOUND*log(arr[L1].rof[d]*RO1);

                for(int i = 0; i < fh->N[d]; i++)
                {
                    switch(d)
                    {
                    case 0:
                        ASSIGN_IND(ind,  i,   j, k);
                        ASSIGN_IND(indR, i+1, j, k);
                        ASSIGN_IND(indL, i-1, j, k);
                        break;
                    case 1:
                        ASSIGN_IND(ind,  k, i,   j);
                        ASSIGN_IND(indR, k, i+1, j);
                        ASSIGN_IND(indL, k, i-1, j);
                        break;
                    case 2:
                        ASSIGN_IND(ind,  j, k, i  );
                        ASSIGN_IND(indR, j, k, i+1);
                        ASSIGN_IND(indL, j, k, i-1);
                        break;
                    }

                    // LOCAL RIEMANN INVARIANTS

                    P1 = SOUND*log(arr[C].ro_fh_n*RO1);
                    P2 = SOUND*log(arr[C].ro_fh*RO1);
                    P3 = SOUND*log(arr[IR].rof[d]*RO1);

                    RCRN = arr[C].u_fh_n[d] + P1;
                    QCRN = arr[C].u_fh_n[d] - P1;

                    RCR  = arr[C].u_fh[d] + P2;
                    QCR  = arr[C].u_fh[d] - P2;

                    RR   = arr[IR].uf[d][d] + P3;
                    QR   = arr[IR].uf[d][d] - P3;

                    VCRN = arr[C].u_fh_n[d1];
                    VCR  = arr[C].u_fh[d1];
                    VR   = arr[IR].uf[d1][d];

                    WCRN = arr[C].u_fh_n[d2];
                    WCR  = arr[C].u_fh[d2];
                    WR   = arr[IR].uf[d2][d];

                    RPR  = RPL;
                    QPR  = QPL;
                    VPR  = VPL;
                    WPR  = WPL;

                    // CHARACTERISTIC SPEEDS
                    L1L = arr[IL].u_fh_n[d] + SOUND;
                    L2L = arr[IL].u_fh_n[d] - SOUND;
                    L3L = arr[IL].u_fh_n[d];
                    L1R = arr[C].u_fh_n[d]  + SOUND;
                    L2R = arr[C].u_fh_n[d]  - SOUND;
                    L3R = arr[C].u_fh_n[d];

                    // EXTRAPOLATION FOR THE LEFT GOING WAVES
                    RPNL = 2.0*RCLN - RPL;
                    QPNL = 2.0*QCLN - QPL;
                    VPNL = 2.0*VCLN - VPL;
                    WPNL = 2.0*WCLN - WPL;

                    // EXTRAPOLATION FOR THE RIGHT GOING WAVES
                    RPNR = 2.0*RCRN - RPR;
                    QPNR = 2.0*QCRN - QPR;
                    VPNR = 2.0*VCRN - VPR;
                    WPNR = 2.0*WCRN - WPR;

                    QQ1 = 2.0*(RCLN - RCL) + DT*L1L*(RPL - RL)/fh->grid.h[C][d];
                    QQ2 = 2.0*(QCLN - QCL) + DT*L2L*(QPL - QL)/fh->grid.h[C][d];
                    QQ3 = 2.0*(VCLN - VCL) + DT*L3L*(VPL - VL)/fh->grid.h[C][d];
                    QQ4 = 2.0*(WCLN - WCL) + DT*L3L*(WPL - WL)/fh->grid.h[C][d];

                    RMAX = DMAX3(RL,RCLN,RPL) + QQ1;    RMIN = DMIN3(RL,RCLN,RPL) + QQ1;
                    QMAX = DMAX3(QL,QCLN,QPL) + QQ2;    QMIN = DMIN3(QL,QCLN,QPL) + QQ2;
                    VMAX = DMAX3(VL,VCLN,VPL) + QQ3;    VMIN = DMIN3(VL,VCLN,VPL) + QQ3;
                    WMAX = DMAX3(WL,WCLN,WPL) + QQ4;    WMIN = DMIN3(WL,WCLN,WPL) + QQ4;

                    if(RPNL > RMAX) RPNL = RMAX;
                    if(RPNL < RMIN) RPNL = RMIN;
                    if(QPNL > QMAX) QPNL = QMAX;
                    if(QPNL < QMIN) QPNL = QMIN;
                    if(VPNL > VMAX) VPNL = VMAX;
                    if(VPNL < VMIN) VPNL = VMIN;
                    if(WPNL > WMAX) WPNL = WMAX;
                    if(WPNL < WMIN) WPNL = WMIN;

                    QQ1 = 2.0*(RCRN - RCR) + DT*L1R*(RR - RPR)/fh->grid.h[C][d];
                    QQ2 = 2.0*(QCRN - QCR) + DT*L2R*(QR - QPR)/fh->grid.h[C][d];
                    QQ3 = 2.0*(VCRN - VCR) + DT*L3R*(VR - VPR)/fh->grid.h[C][d];
                    QQ4 = 2.0*(WCRN - WCR) + DT*L3R*(WR - WPR)/fh->grid.h[C][d];

                    RMAX = DMAX3(RR,RCRN,RPR) + QQ1;    RMIN = DMIN3(RR,RCRN,RPR) + QQ1;
                    QMAX = DMAX3(QR,QCRN,QPR) + QQ2;    QMIN = DMIN3(QR,QCRN,QPR) + QQ2;
                    VMAX = DMAX3(VR,VCRN,VPR) + QQ3;    VMIN = DMIN3(VR,VCRN,VPR) + QQ3;
                    WMAX = DMAX3(WR,WCRN,WPR) + QQ4;    WMIN = DMIN3(WR,WCRN,WPR) + QQ4;

                    if(RPNR > RMAX) RPNR = RMAX;
                    if(RPNR < RMIN) RPNR = RMIN;
                    if(QPNR > QMAX) QPNR = QMAX;
                    if(QPNR < QMIN) QPNR = QMIN;
                    if(VPNR > VMAX) VPNR = VMAX;
                    if(VPNR < VMIN) VPNR = VMIN;
                    if(WPNR > WMAX) WPNR = WMAX;
                    if(WPNR < WMIN) WPNR = WMIN;

                    if((L1L + L1R) >= 0)
                        RPN1 = RPNL;
                    else
                        RPN1 = RPNR;

                    if((L2L + L2R) >= 0)
                        QPN1 = QPNL;
                    else
                        QPN1 = QPNR;

                    if((L3L + L3R) >= 0) {
                        VPN1 = VPNL;
                        WPN1 = WPNL;
                    } else {
                        VPN1 = VPNR;
                        WPN1 = WPNR;
                    }

                    RPN2 = 0.5*(RPNL + RPNR);
                    QPN2 = 0.5*(QPNL + QPNR);
                    VPN2 = 0.5*(VPNL + VPNR);
                    WPN2 = 0.5*(WPNL + WPNR);

                    QQ1 = 0.5*(2.0*(RCRN-RCR+RCLN-RCL) + DT*0.5*(L1L+L1R)*(RR-RPR+RPL-RL)/fh->grid.h[C][d]);
                    QQ2 = 0.5*(2.0*(QCRN-QCR+QCLN-QCL) + DT*0.5*(L2L+L2R)*(QR-QPR+QPL-QL)/fh->grid.h[C][d]);
                    QQ3 = 0.5*(2.0*(VCRN-VCR+VCLN-VCL) + DT*0.5*(L3L+L3R)*(VR-VPR+VPL-VL)/fh->grid.h[C][d]);
                    QQ4 = 0.5*(2.0*(WCRN-WCR+WCLN-WCL) + DT*0.5*(L3L+L3R)*(WR-WPR+WPL-WL)/fh->grid.h[C][d]);

                    RMAX = DMAX3(0.5*(RL+RR), 0.5*(RCLN+RCRN), 0.5*(RPL+RPR)) + QQ1;
                    RMIN = DMIN3(0.5*(RL+RR), 0.5*(RCLN+RCRN), 0.5*(RPL+RPR)) + QQ1;
                    QMAX = DMAX3(0.5*(QL+QR), 0.5*(QCLN+QCRN), 0.5*(QPL+QPR)) + QQ2;
                    QMIN = DMIN3(0.5*(QL+QR), 0.5*(QCLN+QCRN), 0.5*(QPL+QPR)) + QQ2;
                    VMAX = DMAX3(0.5*(VL+VR), 0.5*(VCLN+VCRN), 0.5*(VPL+VPR)) + QQ3;
                    VMIN = DMIN3(0.5*(VL+VR), 0.5*(VCLN+VCRN), 0.5*(VPL+VPR)) + QQ3;
                    WMAX = DMAX3(0.5*(WL+WR), 0.5*(WCLN+WCRN), 0.5*(WPL+WPR)) + QQ4;
                    WMIN = DMIN3(0.5*(WL+WR), 0.5*(WCLN+WCRN), 0.5*(WPL+WPR)) + QQ4;

                    if(RPN2 > RMAX) RPN2 = RMAX;
                    if(RPN2 < RMIN) RPN2 = RMIN;
                    if(QPN2 > QMAX) QPN2 = QMAX;
                    if(QPN2 < QMIN) QPN2 = QMIN;
                    if(VPN2 > VMAX) VPN2 = VMAX;
                    if(VPN2 < VMIN) VPN2 = VMIN;
                    if(WPN2 > WMAX) WPN2 = WMAX;
                    if(WPN2 < WMIN) WPN2 = WMIN;

                    // Blending
                    RPN = blend*RPN1+(1.0 - blend)*RPN2;
                    QPN = blend*QPN1+(1.0 - blend)*QPN2;
                    VPN = blend*VPN1+(1.0 - blend)*VPN2;
                    WPN = blend*WPN1+(1.0 - blend)*WPN2;

                    // New fluxes
                    arr[C].rofn[d]    = fh->FH_dens*exp(0.5*(RPN - QPN)/SOUND);
                    arr[C].ufn[d][d]  = 0.5*(RPN + QPN);
                    arr[C].ufn[d1][d] = VPN;
                    arr[C].ufn[d2][d] = WPN;

                    // Pressure flux
                    switch(fh->eos)
                    {
                    case eos_argon:
                        arr[C].pfn[d] = EOS_D*(1./3.*EOS_A*EOS_A*pow(EOS_E*arr[C].rofn[d] - EOS_B, 3) + EOS_C);
                        break;
                    case eos_spce:
                        arr[C].pfn[d] = arr[C].rofn[d]*(arr[C].rofn[d]*EOS_A + EOS_B) + EOS_C;
                        break;
                    }

                    // For the next step
                    RCLN = RCRN;    QCLN = QCRN;    VCLN = VCRN;    WCLN = WCRN;
                    RCL  = RCR;     QCL  = QCR;     VCL  = VCR;     WCL  = WCR;
                    RL   = RPL;     QL   = QPL;     VL   = VPL;     WL   = WPL;
                    RPL  = RR;      QPL  = QR;      VPL  = VR;      WPL  = WR;
                }
            }
        }
    }
}


void compute_random_stress(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    double FLUCT, TCELL, TRACE;
    double ranv[7];
    double Tcurr, Tf;

    T_INST = 0;

    for(int i = 0; i < fh->Ntot; i++)
    {
        // Fluctuation-Dissipation factor
        FLUCT = sqrt(2.0*FHMD_kB*fh->FH_temp)/sqrt(fh->dt_FH*fh->grid.vol[i]);

        TCELL   = arr[i].ro_fh*(arr[i].u_fh[0]*arr[i].u_fh[0] + arr[i].u_fh[1]*arr[i].u_fh[1] + arr[i].u_fh[2]*arr[i].u_fh[2])*fh->grid.vol[i]/(3.0*FHMD_kB);
        T_INST += TCELL;

        for(int k = 0; k < 7; k++)
            ranv[k] = DRNOR();

        TRACE = (ranv[1] + ranv[2] + ranv[3])/3.0;

        arr[i].rans[0][0] = FLUCT*(sqrt(2.0*MU)*(ranv[1] - TRACE) + sqrt(KAPPA)*ranv[0]);
        arr[i].rans[1][1] = FLUCT*(sqrt(2.0*MU)*(ranv[2] - TRACE) + sqrt(KAPPA)*ranv[0]);
        arr[i].rans[2][2] = FLUCT*(sqrt(2.0*MU)*(ranv[3] - TRACE) + sqrt(KAPPA)*ranv[0]);

        arr[i].rans[0][1] = FLUCT*(sqrt(MU)*ranv[4]);
        arr[i].rans[0][2] = FLUCT*(sqrt(MU)*ranv[5]);
        arr[i].rans[1][2] = FLUCT*(sqrt(MU)*ranv[6]);

        arr[i].rans[1][0] = arr[i].rans[0][1];
        arr[i].rans[2][0] = arr[i].rans[0][2];
        arr[i].rans[2][1] = arr[i].rans[1][2];
    }

    // Compute current temperature and compute dynamic blending factor if enabled
    if(fh->FH_blend < 0)
    {
        blend = 0.0;
        Tcurr = T_INST/(double)(fh->Ntot);

        if(Tcurr > fh->FH_temp)
        {
            Tf = Tcurr/fh->FH_temp;
            blend = 1.0*(Tf - 1.0);
            if(blend > 1.0) blend = 1.0;
        }
    }
}


void swap_var(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    for(int i = 0; i < fh->Ntot; i++)
    {
        for(int d = 0; d < DIM; d++)
        {
            arr[i].rof[d] = arr[i].rofn[d];
            arr[i].pf[d]  = arr[i].pfn[d];

            for(int d1 = 0; d1 < DIM; d1++)
            {
                arr[i].uf[d][d1] = arr[i].ufn[d][d1];
            }
        }
    }
}


void collect_statistics(FHMD *fh)
{
    FH_arrays *arr = fh->arr;

    T_AVG += T_INST/(double)(fh->Ntot);
    NT_AVG++;

    for(int i = 0; i < fh->Ntot; i++)
    {
        RHO_AVG  += arr[i].ro_fh;
        RHO2_AVG += (arr[i].ro_fh - fh->FH_dens)*(arr[i].ro_fh - fh->FH_dens);
        P_AVG    += arr[i].p;
        N_AVG++;

        for(int d = 0; d < DIM; d++)
        {
            U_AVG[d]  += arr[i].u_fh[d];
            U2_AVG[d] += arr[i].u_fh[d]*arr[i].u_fh[d];
        }
    }
}


void update_statistics(FHMD *fh)
{
    const double DN_AVG = (double)(N_AVG);

    avg_rho = RHO_AVG/DN_AVG;
    avg_p   = P_AVG/DN_AVG;
    std_rho = sqrt(fabs(RHO2_AVG/DN_AVG));

    for(int d = 0; d < DIM; d++)
    {
        avg_u[d] = U_AVG[d]/DN_AVG;
        std_u[d] = sqrt(fabs(U2_AVG[d]/DN_AVG - avg_u[d]*avg_u[d]));
        if(std_rho > 0)
            avg_sound[d] = fh->FH_dens*std_u[d]/std_rho;
        else
            avg_sound[d] = 0;
    }

    sound = (avg_sound[0] + avg_sound[1] + avg_sound[2])/3.;
    avg_T = T_AVG/(double)(NT_AVG);
}


void FH_do_single_timestep(FHMD *fh)
{
    compute_random_stress(fh);

    FH_predictor(fh);
    FH_char(fh);
    FH_corrector(fh);

    swap_var(fh);
}


void FH_equilibrate(FHMD *fh)
{
    const double STD_RHO  = sqrt(fh->FH_temp*fh->FH_dens*FHMD_kB/(fh->box_volume/(double)(fh->Ntot)))/SOUND;
    const double STD_U    = sqrt(fh->FH_temp*FHMD_kB/(fh->FH_dens*fh->box_volume/(double)(fh->Ntot)));
    const int    N_OUTPUT = 1000;

    printf(MAKE_BLUE "FHMD: FH equilibration in process...\n\n");
    printf("%8s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %6s\n",
           "Step", "STD rho", "STD Ux", "STD Uy", "STD Uz", "C_s, m/s", "T, K", "<T>, K", "<rho>", "<P>", "<Ux>", "<Uy>", "<Uz>", "blend");
    printf("---------------------------------------------------------------------------------------------------------------------------------------\n");
    printf(MAKE_LIGHT_BLUE "->       %9.4f %9.5f %9.5f %9.5f %9.2f %9.4f %9.4f %9.4f %9.4f      <- Theoretical Values",
           STD_RHO, STD_U, STD_U, STD_U, SOUND*1000.0, fh->FH_temp, fh->FH_temp, fh->FH_dens, P_INIT);
    printf(MAKE_BLUE "\n");

    while(STEP <= fh->FH_equil)
    {
        collect_statistics(fh);

        if(!(STEP % N_OUTPUT))
        {
            update_statistics(fh);

            printf("\r%8d %9.4f %9.5f %9.5f %9.5f %9.2f %9.4f %9.4f %9.4f %9.4f %9.2e %9.2e %9.2e %6.4f",
                   STEP, std_rho, std_u[0], std_u[1], std_u[2], sound*1000.0, T_INST/(double)(fh->Ntot), avg_T,
                   avg_rho, avg_p, avg_u[0], avg_u[1], avg_u[2], blend);

#ifdef FHMD_DEBUG_FH
            printf("\n");
#endif

            fflush(stdout);
        }

        FH_do_single_timestep(fh);

        T += fh->dt_FH;
        STEP++;
    }

    printf("\n---------------------------------------------------------------------------------------------------------------------------------------\n");
    printf(RESET_COLOR "\n");
    fflush(stdout);
}


void define_FH_grid(t_commrec *cr, FHMD *fh)
{
    dvec h0;
    ivec ind;

    for(int d = 0; d < DIM; d++)
        h0[d] = fh->box[d]/(double)(fh->N[d]);

    for(int k = 0; k < NZ; k++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int i = 0; i < NX; i++)
            {
                ASSIGN_IND(ind, i, j, k);

                for(int d = 0; d < DIM; d++)
                {
                    fh->grid.h[C][d] = h0[d];
                    fh->grid.n[C][d] = (double)(ind[d])*h0[d];
                    fh->grid.c[C][d] = fh->grid.n[C][d] + 0.5*h0[d];
                }
            }
        }
    }

    for(int i = 0; i < fh->Ntot; i++)
    {
        fh->grid.vol[i]  = fh->grid.h[i][0]*fh->grid.h[i][1]*fh->grid.h[i][2];
        fh->grid.ivol[i] = 1.0/fh->grid.vol[i];
    }

#ifdef FHMD_DEBUG_GRID
    if(MASTER(cr)) {
        printf(MAKE_YELLOW "FHMD DEBUG: Grid X (i, nodes, centres, steps):\n");
        for(int i = 0; i < fh->N[0]; i++)
            printf("%3d  %8.4f  %8.4f  %8.4f\n", i, fh->grid.n[I3(i,0,0,fh->N)][0], fh->grid.c[I3(i,0,0,fh->N)][0], fh->grid.h[I3(i,0,0,fh->N)][0]);
        printf("%3d  %8.4f     ----      ----\n", fh->N[0], fh->grid.n[I3(fh->N[0],0,0,fh->N)][0]);
        printf(MAKE_YELLOW "FHMD DEBUG: Grid Y (i, nodes, centres, steps):\n");
        for(int i = 0; i < fh->N[1]; i++)
            printf("%3d  %8.4f  %8.4f  %8.4f\n", i, fh->grid.n[I3(0,i,0,fh->N)][1], fh->grid.c[I3(0,i,0,fh->N)][1], fh->grid.h[I3(0,i,0,fh->N)][1]);
        printf("%3d  %8.4f     ----      ----\n", fh->N[1], fh->grid.n[I3(0,fh->N[1],0,fh->N)][1]);
        printf(MAKE_YELLOW "FHMD DEBUG: Grid Z (i, nodes, centres, steps):\n");
        for(int i = 0; i < fh->N[2]; i++)
            printf("%3d  %8.4f  %8.4f  %8.4f\n", i, fh->grid.n[I3(0,0,i,fh->N)][2], fh->grid.c[I3(0,0,i,fh->N)][2], fh->grid.h[I3(0,0,i,fh->N)][2]);
        printf("%3d  %8.4f     ----      ----\n", fh->N[2], fh->grid.n[I3(0,0,fh->N[2],fh->N)][2]);
        printf("==============================================\n");
        printf(RESET_COLOR "\n");
    }
#endif
}
