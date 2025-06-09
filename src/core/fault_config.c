#include "../../include/api/power_system.h"
#include "../../include/core/fault_config.h"

int configure_fault_simulation(NetworkData Net_data_main, AdmittanceMatrix Y_complex_aug_main,
                              double ***Z_AUG_fault, double *fault_start, double *fault_end) {
    
    int fault_enabled = 0;
    int num_buses = Net_data_main.constants[0].LOAD_FLOW_number;
    
    printf("\nSimulate fault? (y/n): ");
    fflush(stdout);
    char fault_ans = 'n';
    scanf(" %c", &fault_ans);

    if (fault_ans == 'y' || fault_ans == 'Y') {
        fault_enabled = 1;

        /* show generator-to-bus mapping */
        int Ng_map = Net_data_main.constants[0].NumberofGens;
        int gen_bus[NG_MAX]; /* use same limit as in solver */
        for (int g = 0; g < Ng_map && g < NG_MAX; ++g) {
            gen_bus[g] = Net_data_main.Load_flow_ps[g].bus_number; /* assumption: order matches */
        }

        printf("Available buses (G = generator bus):\n");
        for (int b = 1; b <= num_buses; ++b) {
            int is_gen = -1;
            for (int g = 0; g < Ng_map; ++g)
                if (gen_bus[g] == b) { is_gen = g; break; }
            if (is_gen >= 0)
                printf("  Bus %d  (Gen %d)\n", b, is_gen);
            else
                printf("  Bus %d\n", b);
        }

        printf("\nAvailable lines (annotated with connected generators):\n");
        int num_lines = Net_data_main.constants[0].Number_of_lines;
        for (int i = 0; i < num_lines; ++i) {
            int b1 = Net_data_main.trans_line_para[i].bus1;
            int b2 = Net_data_main.trans_line_para[i].bus2;
            /* mark if either end is a gen bus */
            int g1 = -1, g2 = -1;
            for (int g = 0; g < Ng_map; ++g) {
                if (gen_bus[g] == b1) g1 = g;
                if (gen_bus[g] == b2) g2 = g;
            }
            char tag1[8]="", tag2[8]="";
            if (g1 >= 0) sprintf(tag1," (G%d)", g1);
            if (g2 >= 0) sprintf(tag2," (G%d)", g2);
            printf("  %2d) Bus %d%s — Bus %d%s\n", i + 1, b1, tag1, b2, tag2);
        }

        /* Ask fault type */
        char fault_type = 'l';
        printf("Fault type — terminal bus (t) or line (l): ");
        fflush(stdout);
        scanf(" %c", &fault_type);

        double fault_duration = 0.0;
        printf("Enter fault start time (s) and duration (s): ");
        fflush(stdout);
        scanf("%lf %lf", fault_start, &fault_duration);
        *fault_end = *fault_start + fault_duration;

        if (fault_type == 't' || fault_type == 'T') {
            int bus_idx = 1;
            printf("Select bus index (1-%d): ", num_buses);
            fflush(stdout);
            scanf("%d", &bus_idx);

            /* Duplicate healthy admittance and ground the selected bus */
            AdmittanceMatrix Y_fault = Y_complex_aug_main; /* shallow copy OK */
            /* Ground via huge shunt admittance to force V ≈ 0 */
            Y_fault.MAT[bus_idx - 1][bus_idx - 1].dat[0] = 0.0;
            Y_fault.MAT[bus_idx - 1][bus_idx - 1].dat[1] = 1e6;

            double **Y_AUG_sp_fault = convertComplexToReal(Y_fault, Net_data_main);
            *Z_AUG_fault = convertComplexToReal(Y_fault, Net_data_main);
            invertMatrix(Y_AUG_sp_fault, 2 * num_buses, *Z_AUG_fault);

            printf("Configured 3-φ bus fault at Bus %d from %.3f to %.3f s\n", 
                   bus_idx, *fault_start, *fault_end);

        } else { /* default to line fault */
            int line_idx = 1;
            printf("Select line index (1-%d): ", num_lines);
            fflush(stdout);
            scanf("%d", &line_idx);

            AdmittanceMatrix Y_fault = createFaultMatrix(line_idx, Net_data_main);
            double **Y_AUG_sp_fault = convertComplexToReal(Y_fault, Net_data_main);
            *Z_AUG_fault = convertComplexToReal(Y_fault, Net_data_main);
            invertMatrix(Y_AUG_sp_fault, 2 * num_buses, *Z_AUG_fault);

            int b1 = Net_data_main.trans_line_para[line_idx - 1].bus1;
            int b2 = Net_data_main.trans_line_para[line_idx - 1].bus2;
            printf("Configured line fault %d (Bus %d — Bus %d) from %.3f to %.3f s\n",
                   line_idx, b1, b2, *fault_start, *fault_end);
        }
    }
    
    return fault_enabled;
} 