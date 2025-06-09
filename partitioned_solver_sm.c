#include <stdio.h>
#include <math.h>
#include "read_header.h"
#include "data_types.h"

/*
 * External functions (already implemented in func2.c) that we reuse here.
 * We only declare them so the compiler knows their signatures.
 */

double Power_elc(double Ed2,double id,double Eq2,double iq,double Xd2,double Xq2);
void Network_solver(    COMB_STRUCT DATA,
                        int Number_of_generators,
                        int Number_of_buses,
                        NW_STATES * pointer_nw_vector,
                        STATES * ptr_X_vector,
                        double ** Z_AUG,
                        double * ptr_iq,
                        double * ptr_id,
                        double * ptr_Vt);

void current_update_new(STATES X_iter,
                        double vd_iter,
                        double vq_iter,
                        double * ptr_id,
                        double * ptr_iq,
                        double Xd2_iter);

/* -------------------------------------------------------------------------
 *  Internal generator-state derivative and Heun stepper implementation     *
 * -------------------------------------------------------------------------*/

#ifndef OMEGA_BASE
#define OMEGA_BASE (2.0 * M_PI * 60.0)   /* rad/s at 60-Hz base */
#endif

typedef struct {
    double dEq1, dEd1, dEq2, dEd2;
    double ddelta, dslip;
    double dEfd;
} GEN_DERIV;

static GEN_DERIV compute_deriv(int g,
                               const COMB_STRUCT *D,
                               const STATES *X,
                               double id,
                               double iq,
                               double Pm,
                               double Pe,
                               double D_by_M,
                               double prop_err)
{
    const MCD  *Gp = &D->Generator_ps[g];
    const EXCT *Ep = &D->exciter_ps[g];

    double tdo1 = Gp->Tdo1, tqo1 = Gp->Tqo1;
    double tdo2 = Gp->Tdo2, tqo2 = Gp->Tqo2;

    double xd0 = Gp->Xd0, xq0 = Gp->Xq0;
    double xd1 = Gp->Xd1, xq1 = Gp->Xq1;
    double xd2 = Gp->Xd2, xq2 = Gp->Xq2;

    double one_by_2H = 1.0 / (2.0 * Gp->H);

    GEN_DERIV f;

    /* flux-linkage dynamics */
    f.dEq1 = (X->Efd - X->Eq_das + id * (xd0 - xd1)) / tdo1;
    f.dEd1 = -(X->Ed_das + iq * (xq0 - xq1)) / tqo1;
    f.dEq2 = (X->Eq_das - X->Eq_das_das + id * (xd1 - xd2)) / tdo2;
    f.dEd2 = (X->Ed_das - X->Ed_das_das + iq * (xq2 - xq1)) / tqo2;

    /* mechanical ‑ electrical swing */
    double torque_imb = Pm - Pe;
    f.dslip  = torque_imb * one_by_2H - D_by_M * X->slip;
    f.ddelta = OMEGA_BASE * X->slip;

    /* AVR first-order model */
    f.dEfd   = (Ep->ka * prop_err - X->Efd) / Ep->ta;

    return f;
}

static void advance_state(STATES *S, const GEN_DERIV k, double h) {
    S->Eq_das       += h * k.dEq1;
    S->Ed_das       += h * k.dEd1;
    S->Eq_das_das   += h * k.dEq2;
    S->Ed_das_das   += h * k.dEd2;
    S->delta        += h * k.ddelta;
    S->slip         += h * k.dslip;
    S->Efd          += h * k.dEfd;
}

static void rk4_step(int g,
                     const COMB_STRUCT *D,
                     STATES *X,
                     double dt,
                     double id,
                     double iq,
                     double Pm,
                     double Pe,
                     double D_by_M,
                     double prop_err)
{
    GEN_DERIV k1 = compute_deriv(g, D, X, id, iq, Pm, Pe, D_by_M, prop_err);

    STATES X2 = *X;
    advance_state(&X2, k1, 0.5*dt);
    GEN_DERIV k2 = compute_deriv(g, D, &X2, id, iq, Pm, Pe, D_by_M, prop_err);

    STATES X3 = *X;
    advance_state(&X3, k2, 0.5*dt);
    GEN_DERIV k3 = compute_deriv(g, D, &X3, id, iq, Pm, Pe, D_by_M, prop_err);

    STATES X4 = *X;
    advance_state(&X4, k3, dt);
    GEN_DERIV k4 = compute_deriv(g, D, &X4, id, iq, Pm, Pe, D_by_M, prop_err);

    /* weighted sum */
    X->Eq_das       += dt/6.0 * (k1.dEq1 + 2*k2.dEq1 + 2*k3.dEq1 + k4.dEq1);
    X->Ed_das       += dt/6.0 * (k1.dEd1 + 2*k2.dEd1 + 2*k3.dEd1 + k4.dEd1);
    X->Eq_das_das   += dt/6.0 * (k1.dEq2 + 2*k2.dEq2 + 2*k3.dEq2 + k4.dEq2);
    X->Ed_das_das   += dt/6.0 * (k1.dEd2 + 2*k2.dEd2 + 2*k3.dEd2 + k4.dEd2);
    X->delta        += dt/6.0 * (k1.ddelta + 2*k2.ddelta + 2*k3.ddelta + k4.ddelta);
    X->slip         += dt/6.0 * (k1.dslip  + 2*k2.dslip  + 2*k3.dslip  + k4.dslip);
    X->Efd          += dt/6.0 * (k1.dEfd   + 2*k2.dEfd   + 2*k3.dEfd   + k4.dEfd);

    /* clamp Efd and leave delta unwrapped as before */
}

/*****************************************************************************************
 *  Simplified partitioned time-domain solver with internal Heun integration             *
 *                                                                                       *
 *  – Works for any number of generators / buses provided in All_data                    *
 *  – Uses a single, healthy impedance matrix (Z_AUG_healthy) for the whole run          *
 *  – No fault logic; no copy-back predictor; no pointer gymnastics                      *
 *  – Outputs only per-generator CSV files in sim/ directory (no main CSV)              *
 *                                                                                       *
 *  Author: chatGPT refactor                                                             *
 *****************************************************************************************/
#define NG_MAX 20
static double xi[NG_MAX] = {0};   /* integrator */

void partitioned_solver_sm(
        COMB_STRUCT All_data,
        INITIAL   *Initial_state,
        double     del_t,
        double     END_TIME,
        double   **Z_AUG_healthy,
        double   **Z_AUG_fault,
        int        fault_enabled,
        double     fault_start,
        double     fault_end,
        int        TARGET_G,
        double     VREF_STEP_TIME,
        double     VREF_STEP_DELTA,
        double     PM_STEP_TIME,
        double     PM_STEP_DELTA)
{
    const int Ng = All_data.constants[0].NumberofGens;
    const int Nb = All_data.constants[0].LOAD_FLOW_number;

    /* ----------------------- local working buffers ---------------------------------- */
    STATES       X[Ng];            /* generator state vector                        */
    NW_STATES    NW[Nb];           /* network bus voltages (dq)                     */
    double       id[Ng], iq[Ng];   /* stator currents                               */
    double       Vt[Ng];
    double       vq_g[Ng], vd_g[Ng];
    double       Vref[Ng];         /* AVR reference voltage                         */
    double       Pe[Ng];           /* electrical power                              */
    double       Pm[Ng];           /* mechanical power (dynamic with governor)    */
    double       Pref[Ng];         /* speed reference for governor (initial Pm0)  */
    double       D_by_M[Ng];       /* damping per unit inertia                      */

    /* ----------------------- initialisation ---------------------------------------- */
    for (int g = 0; g < Ng; ++g) {
        /* mechanical parameters */
        D_by_M[g] = 2.0 + 0.5 * g;   /* Much higher damping */

        /* copy initial state */
        X[g].Ed_das      = Initial_state[g].Ed_dash_0.dat[0];
        X[g].Eq_das      = Initial_state[g].Eq_dash_0.dat[0];
        X[g].Ed_das_das  = Initial_state[g].Ed_das_das;
        X[g].Eq_das_das  = Initial_state[g].Eq_das_das;
        X[g].delta       = Initial_state[g].delta_0;
        X[g].slip        = 0.0;
        X[g].Efd         = Initial_state[g].Efd_0.dat[0];
        X[g].E_dummy     = 0.0;   /* will be made consistent after first curr. calc */

        /* mechanical/electrical power */
        Pm[g]   = Initial_state[g].Pm_0.dat[0];
        Pref[g] = Pm[g];   /* store as reference set-point */
        Pe[g] = Initial_state[g].Pe_0.dat[0];

        /* reference voltage calculated from initial conditions */
        double Vq0 = Initial_state[g].vq_0.dat[0];
        double Vd0 = Initial_state[g].vd_0.dat[0];
        Vt[g]      = sqrt(Vq0*Vq0 + Vd0*Vd0);
        Vref[g]    = Vt[g] + (Initial_state[g].Efd_0.dat[0] / All_data.exciter_ps[g].ka);

        /* placeholder currents (will be overwritten after first network solve) */
        iq[g] = Initial_state[g].iq_0.dat[0];
        id[g] = Initial_state[g].id_0.dat[0];
    }

    /* after first Network_solver(), before the main RK loop */
    for (int g = 0; g < Ng; ++g) {
        Pm[g]   = Pe[g];   /* EXACT match */
        Pref[g] = Pm[g];
    }

    /* open a CSV file per generator */
    FILE *gen_csv[Ng];
    for (int g = 0; g < Ng; ++g) {
        char fname[64];
        sprintf(fname, "sim/gen%d.csv", g);
        gen_csv[g] = fopen(fname, "w");
        if (!gen_csv[g]) { perror("open gen csv"); exit(1); }
        fprintf(gen_csv[g],
                "time,delta,slip,Efd,Eq2,Ed2,Eq1,Ed1,E_dummy,Vt,Vref,vq,vd,id_0,iq_0,mech_power,elec_power\n");
    }

    /* ----------------------- main time-march --------------------------------------- */
    for (double t = 0.0; t < END_TIME; t += del_t) {

        /* 1. network solution with latest states */
        double **Z_use = (fault_enabled && t >= fault_start && t < fault_end) ? Z_AUG_fault : Z_AUG_healthy;

        Network_solver(All_data, Ng, Nb,
                       NW,               /* output dq bus voltages          */
                       X,                /* generator states (uses delta)   */
                       Z_use,            /* impedance matrix (healthy/fault)*/
                       iq, id,           /* currents arrays (in/out)        */
                       Vt);              /* terminal voltage magnitude out  */

        /* 2. update currents based on the fresh bus voltages */
        for (int g = 0; g < Ng; ++g) {
            double vq = NW[g].VQ * cos(X[g].delta) + NW[g].VD * sin(X[g].delta);
            double vd = NW[g].VD * cos(X[g].delta) - NW[g].VQ * sin(X[g].delta);
            vq_g[g]=vq; vd_g[g]=vd;
            current_update_new(X[g], vd, vq, &id[g], &iq[g], All_data.Generator_ps[g].Xd2);
            /* keep E_dummy consistent */
            X[g].E_dummy = -(All_data.Generator_ps[g].Xq2 - All_data.Generator_ps[g].Xd2) * iq[g];
        }

        /* 3. electrical power */
        for (int g = 0; g < Ng; ++g) {
            Pe[g] = Power_elc(X[g].Ed_das_das, id[g],
                               X[g].Eq_das_das, iq[g],
                               All_data.Generator_ps[g].Xd2,
                               All_data.Generator_ps[g].Xq2);
        }

        /* --------------------------------------------------------------
         *  First-step initialisation tweak: force Pm = Pe so that the
         *  accelerating torque (Pm−Pe) starts at zero.  This removes the
         *  large initial swing in slip/δ that came from a mismatch
         *  between load-flow electrical power and assumed mechanical
         *  set-point.
         * --------------------------------------------------------------*/
        static int first_iter_done = 0;
        if (!first_iter_done) {
            for (int g = 0; g < Ng; ++g) {
                Pm[g]   = Pe[g];
                Pref[g] = Pm[g];
            }
            first_iter_done = 1;
        }

        /* 3b. simple governor: first-order droop model  ---------------------------- */
        const double R_GOV   = 0.01;   /* 1% droop (tighter) */
        const double TGOV    = 0.2;    /* faster response */
        const double KI_GOV  = 1.0;    /* stronger integral */

        for (int g = 0; g < Ng; ++g) {
            double speed_err = -X[g].slip;      /* positive if machine slow */

            /* PI droop: Pref + Kp*err + integral */
            xi[g] += KI_GOV * speed_err * del_t;
            double Pm_cmd = Pref[g] + speed_err / R_GOV + xi[g];

            /* first-order actuator lag */
            Pm[g] += (Pm_cmd - Pm[g]) * (del_t / TGOV);

            if (Pm[g] < 0.0)  Pm[g] = 0.0;
            if (Pm[g] > 1.5) Pm[g] = 1.5;
        }

        /* 4. integrate generator dynamics (internal Heun step) */
        for (int g = 0; g < Ng; ++g) {
            double prop_err = Vref[g] - Vt[g];
            rk4_step(g, &All_data, &X[g], del_t,
                      id[g], iq[g],
                      Pm[g], Pe[g],
                      D_by_M[g], prop_err);
        }

        /* 5. dump to per-generator CSV files only */
        for (int g = 0; g < Ng; ++g) {
            /* wrap to ±10π instead of ±π */
            double delta_out = fmod(X[g].delta + 10.0*M_PI, 20.0*M_PI) - 10.0*M_PI;
            fprintf(gen_csv[g], "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                    t,
                    delta_out,
                    X[g].slip,
                    X[g].Efd,
                    X[g].Eq_das_das,
                    X[g].Ed_das_das,
                    X[g].Eq_das,
                    X[g].Ed_das,
                    X[g].E_dummy,
                    Vt[g], Vref[g], vq_g[g], vd_g[g],
                    id[g], iq[g],
                    Pm[g], Pe[g]);
        }

        /* 6. step changes for selected generator */
        static int vref_done = 0, pm_done = 0;
        if (!vref_done && t >= VREF_STEP_TIME) {
            Vref[TARGET_G] += VREF_STEP_DELTA;
            printf("[solver] Gen %d Vref step %+0.4f pu at t=%.1f s\n",
                   TARGET_G, VREF_STEP_DELTA, t);
            vref_done = 1;
        }

        if (!pm_done && t >= PM_STEP_TIME) {
            Pm[TARGET_G] += PM_STEP_DELTA;
            printf("[solver] Gen %d Pm step %+0.4f pu at t=%.1f s\n",
                   TARGET_G, PM_STEP_DELTA, t);
            pm_done = 1;
        }
    }

    /* close per-gen files */
    for (int g = 0; g < Ng; ++g) fclose(gen_csv[g]);
} 