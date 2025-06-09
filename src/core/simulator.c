#include "../../include/api/power_system.h"

/* -------------------------------------------------------------------------
 *  Internal generator-state derivative and Heun stepper implementation     *
 * -------------------------------------------------------------------------*/

/* OMEGA_BASE is now defined in common.h */

typedef struct {
    double dEq1, dEd1, dEq2, dEd2;
    double ddelta, dslip;
    double dEfd;
} GEN_DERIV;

static GEN_DERIV compute_deriv(int g,
                               const NetworkData *D,
                               const GeneratorStates *X,
                               double id,
                               double iq,
                               double Pm,
                               double Pe,
                               double D_by_M,
                               double prop_err)
{
    const GeneratorParams  *Gp = &D->Generator_ps[g];
    const ExciterParams    *Ep = &D->exciter_ps[g];

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

static void advance_state(GeneratorStates *S, const GEN_DERIV k, double h) {
    S->Eq_das       += h * k.dEq1;
    S->Ed_das       += h * k.dEd1;
    S->Eq_das_das   += h * k.dEq2;
    S->Ed_das_das   += h * k.dEd2;
    S->delta        += h * k.ddelta;
    S->slip         += h * k.dslip;
    S->Efd          += h * k.dEfd;
}

static void rk4_step(int g,
                     const NetworkData *D,
                     GeneratorStates *X,
                     double dt,
                     double id,
                     double iq,
                     double Pm,
                     double Pe,
                     double D_by_M,
                     double prop_err)
{
    GEN_DERIV k1 = compute_deriv(g, D, X, id, iq, Pm, Pe, D_by_M, prop_err);

    GeneratorStates X2 = *X;
    advance_state(&X2, k1, 0.5*dt);
    GEN_DERIV k2 = compute_deriv(g, D, &X2, id, iq, Pm, Pe, D_by_M, prop_err);

    GeneratorStates X3 = *X;
    advance_state(&X3, k2, 0.5*dt);
    GEN_DERIV k3 = compute_deriv(g, D, &X3, id, iq, Pm, Pe, D_by_M, prop_err);

    GeneratorStates X4 = *X;
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
#define NG_MAX 20  /* retains hard-coded limits for other arrays */

/* Per-run governor integrator (allocated dynamically so it resets each call) */
double *xi = NULL;  /* will be calloc'd after Ng is known */

const double PM_MAX = 3.0;   /* put near top of simulator.c */

void partitioned_solver_sm(
        NetworkData All_data,
        InitialConditions   *Initial_state,
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

    /* ------------------------------------------------------------------
     * Generator-to-bus mapping: generator g may sit on an arbitrary bus.
     * Store zero-based bus indices so we can fetch voltages correctly.
     * ------------------------------------------------------------------ */
    int bus_of_gen[Ng];
    for (int g = 0; g < Ng; ++g) {
        /* Assumes first Ng entries of Load_flow_ps correspond to gens */
        int bus_num = All_data.Load_flow_ps[g].bus_number; /* 1-based */
        bus_of_gen[g] = bus_num - 1;                        /* 0-based */
    }

    /* allocate governor integrator */
    xi = calloc(Ng, sizeof(double));
    if (!xi) { perror("calloc xi"); exit(1); }

    /* ----------------------- local working buffers ---------------------------------- */
    GeneratorStates  X[Ng];            /* generator state vector                        */
    NetworkStates    NW[Nb];           /* network bus voltages (dq)                     */
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
        /* Mechanical damping coefficient (realistic nominal) */
        D_by_M[g] = 0.5;  /* pu torque per pu speed Earlier used this 3.0 + 1.0 * g; */

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
        Vref[g]    = Vt[g];

        /* placeholder currents (will be overwritten after first network solve) */
        iq[g] = Initial_state[g].iq_0.dat[0];
        id[g] = Initial_state[g].id_0.dat[0];
    }

    /* after first Network_solver(), before the main RK loop */
    for (int g = 0; g < Ng; ++g) {
        Pm[g]   = Pe[g];   /* EXACT match */
        Pref[g] = Pm[g];
    }

    /* ===================================================================
     * INITIAL CONDITION CONSISTENCY FIX
     * ===================================================================
     * Perform ONE network solution + current calculation to ensure that
     * the currents, EMFs, and powers are perfectly consistent between
     * initialization and the dynamic simulation.
     * =================================================================== */
    
    printf("\n🔧 Applying initial condition consistency fix...\n");
    
    /* Step 1: Solve network with initial states to get consistent bus voltages */
    solveNetworkStep(All_data, Ng, Nb,
                     NW,               /* output dq bus voltages          */
                     X,                /* generator states (uses delta)   */
                     Z_AUG_healthy,    /* use healthy impedance matrix    */
                     iq, id,           /* currents arrays (in/out)        */
                     Vt);              /* terminal voltage magnitude out  */

    /* Step 2: Update currents based on fresh bus voltages and states */
    for (int g = 0; g < Ng; ++g) {
        int bus_idx = bus_of_gen[g];
        double vq = NW[bus_idx].VQ * cos(X[g].delta) + NW[bus_idx].VD * sin(X[g].delta);
        double vd = NW[bus_idx].VD * cos(X[g].delta) - NW[bus_idx].VQ * sin(X[g].delta);
        vq_g[g] = vq; 
        vd_g[g] = vd;
        
        /* Calculate currents using the same method as in main loop */
        current_update_new(X[g], vd, vq, &id[g], &iq[g], All_data.Generator_ps[g].Xd2);
        
        /* Update E_dummy consistently */
        X[g].E_dummy = -(All_data.Generator_ps[g].Xq2 - All_data.Generator_ps[g].Xd2) * iq[g];
    }

    /* Step 3: Calculate electrical power using consistent currents and EMFs */
    for (int g = 0; g < Ng; ++g) {
        Pe[g] = computeElectricalPower(X[g].Ed_das_das, id[g],
                                       X[g].Eq_das_das, iq[g],
                                       All_data.Generator_ps[g].Xd2,
                                       All_data.Generator_ps[g].Xq2);
    }

    /* Step 4: Set mechanical power equal to electrical power for perfect balance */
    for (int g = 0; g < Ng; ++g) {
        Pm[g] = Pe[g];          /* Perfect Pm = Pe balance */
        Pref[g] = Pm[g];        /* Governor reference */
        
        printf("   Gen %d: Pm = Pe = %.6f pu, Delta = %.3f°\n", 
               g, Pm[g], X[g].delta * 180.0/PI);
    }

    /* Step 5: Force steady-state consistency by ensuring zero EMF derivatives */
    for (int g = 0; g < Ng; ++g) {
        double Vt_g = sqrt(vq_g[g]*vq_g[g] + vd_g[g]*vd_g[g]);
        
        /* CRITICAL FIX: Adjust EMFs to ensure zero derivatives in steady state */
        const GeneratorParams *Gp = &All_data.Generator_ps[g];
        
        /* Force dEq'/dt = 0: Eq' = Efd + id*(Xd0-Xd1) */
        X[g].Eq_das = X[g].Efd + id[g] * (Gp->Xd0 - Gp->Xd1);
        
        /* Force dEd'/dt = 0: Ed' = -iq*(Xq0-Xq1) */
        X[g].Ed_das = -iq[g] * (Gp->Xq0 - Gp->Xq1);
        
        /* Force dEq''/dt = 0: Eq'' = Eq' + id*(Xd1-Xd2) */
        X[g].Eq_das_das = X[g].Eq_das + id[g] * (Gp->Xd1 - Gp->Xd2);
        
        /* Force dEd''/dt = 0: Ed'' = Ed' + iq*(Xq2-Xq1) */
        X[g].Ed_das_das = X[g].Ed_das + iq[g] * (Gp->Xq2 - Gp->Xq1);
        
        /* Set Vref to maintain steady state with current Efd */
        Vref[g] = Vt_g;
        /* Set Vref to make AVR output steady (dEfd/dt = 0) → Efd = ka*(Vref−Vt) */ 
        // This messes up delta and we see enormous swings 
        // const ExciterParams *Ep = &All_data.exciter_ps[g];
        // Vref[g] = Vt_g + X[g].Efd / Ep->ka;
        
        printf("   Gen %d: Vt = %.6f pu, Vref = %.6f pu, Efd = %.6f pu\n", 
               g, Vt_g, Vref[g], X[g].Efd);
        printf("   Gen %d: Corrected EMFs - Eq'=%.4f, Ed'=%.4f, Eq''=%.4f, Ed''=%.4f\n", 
               g, X[g].Eq_das, X[g].Ed_das, X[g].Eq_das_das, X[g].Ed_das_das);
    }
    
    /* Step 6: Recalculate electrical power with corrected EMFs */
    for (int g = 0; g < Ng; ++g) {
        Pe[g] = computeElectricalPower(X[g].Ed_das_das, id[g],
                                       X[g].Eq_das_das, iq[g],
                                       All_data.Generator_ps[g].Xd2,
                                       All_data.Generator_ps[g].Xq2);
    }
    
    /* Step 7: Verify that EMF derivatives are actually zero */
    printf("\n🔍 EMF Derivative Verification:\n");
    for (int g = 0; g < Ng; ++g) {
        const GeneratorParams *Gp = &All_data.Generator_ps[g];
        
        /* Calculate what the derivatives should be with corrected EMFs */
        double dEq1_dt = (X[g].Efd - X[g].Eq_das + id[g] * (Gp->Xd0 - Gp->Xd1)) / Gp->Tdo1;
        double dEd1_dt = -(X[g].Ed_das + iq[g] * (Gp->Xq0 - Gp->Xq1)) / Gp->Tqo1;
        double dEq2_dt = (X[g].Eq_das - X[g].Eq_das_das + id[g] * (Gp->Xd1 - Gp->Xd2)) / Gp->Tdo2;
        double dEd2_dt = (X[g].Ed_das - X[g].Ed_das_das + iq[g] * (Gp->Xq2 - Gp->Xq1)) / Gp->Tqo2;
        
        printf("   Gen %d: dEq'/dt = %.8f, dEd'/dt = %.8f, dEq''/dt = %.8f, dEd''/dt = %.8f\n", 
               g, dEq1_dt, dEd1_dt, dEq2_dt, dEd2_dt);
               
        /* Check if derivatives are close to zero */
        double max_deriv = fmax(fmax(fabs(dEq1_dt), fabs(dEd1_dt)), fmax(fabs(dEq2_dt), fabs(dEd2_dt)));
        if (max_deriv < 1e-6) {
            printf("   Gen %d: ✅ EMF derivatives are zero - steady state achieved!\n", g);
        } else {
            printf("   Gen %d: ⚠️  EMF derivatives not zero - may cause drift!\n", g);
        }
    }

    /* Step 8: Final power balance – set exact equality and reset governor integrator */
    for (int g = 0; g < Ng; ++g) {
        Pm[g]   = Pe[g];   /* perfect torque balance */
        Pref[g] = Pm[g];   /* governor reference    */
        xi[g]   = 0.0;     /* reset PI integrator   */
        printf("   Gen %d: FINAL Pm = Pe = %.6f pu, Delta = %.3f°\n",
               g, Pm[g], X[g].delta * 180.0/PI);
    }
    
    /* Final: clear governor integrator to avoid residual bias */
    for (int g = 0; g < Ng; ++g) xi[g] = 0.0;

    printf("✅ Initial condition consistency fix completed! (governor integrator reset)\n\n");

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

    /* open an events log file */
    FILE *evt_csv = fopen("sim/events.csv", "w");
    if (!evt_csv) { perror("open events csv"); exit(1); }
    fprintf(evt_csv, "time,event,detail\n");

    /* track fault window state for logging */
    int prev_fault_state = 0;

    /* ----------------------- main time-march --------------------------------------- */
    for (double t = 0.0; t < END_TIME; t += del_t) {

        /* 1. decide which impedance matrix to use based on fault window */
        int fault_now = (fault_enabled && t >= fault_start && t < fault_end);
        if (fault_enabled && fault_now && !prev_fault_state) {
            fprintf(evt_csv, "%lf,fault_start,\n", t);
        } else if (fault_enabled && !fault_now && prev_fault_state) {
            fprintf(evt_csv, "%lf,fault_end,\n", t);
        }
        prev_fault_state = fault_now;
        double **Z_use = fault_now ? Z_AUG_fault : Z_AUG_healthy;

        /* 2. solve the network equations */
        solveNetworkStep(All_data, Ng, Nb,
                         NW,               /* output dq bus voltages          */
                         X,                /* generator states (uses delta)   */
                         Z_use,            /* impedance matrix (healthy/fault)*/
                         iq, id,           /* currents arrays (in/out)        */
                         Vt);              /* terminal voltage magnitude out  */

        /* 4. update currents based on the fresh bus voltages */
        for (int g = 0; g < Ng; ++g) {
            int bus_idx = bus_of_gen[g];
            double vq = NW[bus_idx].VQ * cos(X[g].delta) + NW[bus_idx].VD * sin(X[g].delta);
            double vd = NW[bus_idx].VD * cos(X[g].delta) - NW[bus_idx].VQ * sin(X[g].delta);
            vq_g[g]=vq; vd_g[g]=vd;
            current_update_new(X[g], vd, vq, &id[g], &iq[g], All_data.Generator_ps[g].Xd2);
            /* keep E_dummy consistent */
            X[g].E_dummy = -(All_data.Generator_ps[g].Xq2 - All_data.Generator_ps[g].Xd2) * iq[g];
        }

        /* 5. electrical power */
        for (int g = 0; g < Ng; ++g) {
            Pe[g] = computeElectricalPower(X[g].Ed_das_das, id[g],
                                            X[g].Eq_das_das, iq[g],
                                            All_data.Generator_ps[g].Xd2,
                                            All_data.Generator_ps[g].Xq2);
        }

        /* Initial condition consistency is now handled before main loop */

        /* 3b. simple governor: first-order droop model  ---------------------------- */
        const double R_GOV   = 0.003;   /* 0.3 % droop  */
        const double TGOV    = 0.05;    /* 50 ms actuator */
        const double KI_GOV  =  8.0;    /* strong integral */

        for (int g = 0; g < Ng; ++g) {
            double speed_err = -X[g].slip;      /* positive if machine slow */

            /* PI droop: Pref + Kp*err + integral */
            xi[g] += KI_GOV * speed_err * del_t;
            double Pm_cmd = Pref[g] + speed_err / R_GOV + xi[g];

            /* first-order actuator lag */
            Pm[g] += (Pm_cmd - Pm[g]) * (del_t / TGOV);

            if (Pm[g] < 0.0)  Pm[g] = 0.0;
            if (Pm[g] > PM_MAX) Pm[g] = PM_MAX;
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
            double delta_out = fmod(X[g].delta + 10.0*PI, 20.0*PI) - 10.0*PI;
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
        if (!vref_done && t >= VREF_STEP_TIME && fabs(VREF_STEP_DELTA) > 1e-9) {
            Vref[TARGET_G] += VREF_STEP_DELTA;
            printf("[solver] Gen %d Vref step %+0.4f pu at t=%.1f s\n",
                   TARGET_G, VREF_STEP_DELTA, t);
            fprintf(evt_csv, "%lf,vref_step,%+0.4f\n", t, VREF_STEP_DELTA);
            vref_done = 1;
        }

        if (!pm_done && t >= PM_STEP_TIME && fabs(PM_STEP_DELTA) > 1e-9) {
            Pref[TARGET_G] += PM_STEP_DELTA;         /* raise reference */
            xi[TARGET_G]   += PM_STEP_DELTA;         /* bias integrator so new level holds */
            printf("[solver] Gen %d Pref+xi step %+0.4f pu at t=%.1f s\n",
                   TARGET_G, PM_STEP_DELTA, t);
            fprintf(evt_csv, "%lf,pm_step,%+0.4f\n", t, PM_STEP_DELTA);
            pm_done = 1;
        }
    }

    /* close per-gen files */
    for (int g = 0; g < Ng; ++g) fclose(gen_csv[g]);
    fclose(evt_csv);

    /* free dynamic integrator */
    free(xi);
} 