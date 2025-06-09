#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../include/io/user_input.h"
#include "../../include/common/data_types.h"

void get_simulation_parameters(int *target_g, double *vref_step_delta, 
                              double *pm_step_delta, double *vref_step_time, 
                              double *pm_step_time, InitialConditions *Initial_state_main) {
    
    /* Set defaults */
    *vref_step_time = 200.0;  /* default 100 s */
    *vref_step_delta = 0.05;   /* default +0.05 pu */
    *pm_step_time = 200.0;     /* default same as Vref */
    *pm_step_delta = 0.0;      /* default 0 (no change) */
    *target_g = 0;             /* default generator index */

    /* ------------------------------------------------------------------
     * Flexible user-input parser
     * ------------------------------------------------------------------
     * Accept up to three whitespace-separated numbers in ANY of the forms:
     *   1) "<ΔVref>"                       (Pm step = 0, gen = 0)
     *   2) "<ΔVref> <ΔPm>"                 (gen = 0)
     *   3) "<gen> <ΔVref>"                 (Pm step = 0)
     *   4) "<gen> <ΔVref> <ΔPm>"           (fully specified)
     *   5) empty line                       (all defaults)
     * The first token is treated as generator index IFF it has no decimal
     * point; otherwise it is interpreted as ΔVref.
     * ------------------------------------------------------------------ */

    char line[128] = {0};
    printf("Enter generator index, Vref Δ, Pm Δ  (⏎ for defaults): ");
    fflush(stdout);
    if (!fgets(line, sizeof(line), stdin)) {
        line[0] = '\0'; /* treat as empty */
    }

    int   gen_local    = 0;
    double vref_local  = 0.0;
    double pm_local    = 0.0;

    /* tokenise */
    char *tok = strtok(line, " \t\n");
    char *tokens[3];
    int ntok = 0;
    while (tok && ntok < 3) {
        tokens[ntok++] = tok;
        tok = strtok(NULL, " \t\n");
    }

    if (ntok == 0) {
        /* keep defaults */
    } else {
        /* helper to detect integer */
        int treat_first_as_gen = 0;
        if (strchr(tokens[0], '.') == NULL && strchr(tokens[0], 'e') == NULL && strchr(tokens[0], 'E') == NULL) {
            treat_first_as_gen = 1;
        }

        int idx = 0;
        if (treat_first_as_gen) {
            gen_local = atoi(tokens[idx++]);
        }

        if (idx < ntok) vref_local = atof(tokens[idx++]);
        if (idx < ntok) pm_local   = atof(tokens[idx++]);
    }

    *target_g         = gen_local;
    *vref_step_delta  = vref_local;
    *pm_step_delta    = pm_local;

    /* Step timing */
    printf("Enter step time in seconds (⏎ for default 100): ");
    fflush(stdout);
    if (!fgets(line, sizeof(line), stdin) || line[0] == '\n') {
        *vref_step_time = 100.0;
    } else {
        *vref_step_time = atof(line);
        if (*vref_step_time <= 0.0) *vref_step_time = 100.0;
    }

    *pm_step_time = *vref_step_time; /* same moment for both steps */

    /* Show current reference values */
    printf("Selected Generator %d | Initial Vt ≈ %.4f pu | Initial Pm ≈ %.4f pu\n",
           *target_g,
           sqrt(pow(Initial_state_main[*target_g].vq_0.dat[0],2)+pow(Initial_state_main[*target_g].vd_0.dat[0],2)),
           Initial_state_main[*target_g].Pm_0.dat[0]);
    printf("Scheduled: Vref Δ = %+0.4f pu, Pm Δ = %+0.4f pu at t = %.1f s\n",
           *vref_step_delta, *pm_step_delta, *vref_step_time);

    /* Inform user about governor mechanical-power limits so they understand
       why an excessively large step may be clipped. Keep the hard-coded
       ceiling value consistent with simulator.c (currently 1.5 pu). */
    printf("[info] Governor model clamps mechanical power to the range 0.0 – 1.5 pu.\n");
    printf("       Steps that would push Pm outside this range will appear as a spike then\n");
    printf("       settle at the nearest limit. Adjust ΔPm accordingly for a lasting step.\n");
} 