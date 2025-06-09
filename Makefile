.PHONY: all clean test analysis plots_gen

# Default: build the executable
all: build

# Compile everything
build:
	@echo "→ Building simulation executable"
	@mkdir -p report_plot sim
	@gcc -g -Wall -c main.c functions.c read_fn.c invert.c func2.c partitioned_solver_sm.c
	@gcc main.o functions.o read_fn.o invert.o func2.o partitioned_solver_sm.o -lgsl -lgslcblas -lumfpack -lm -o test
	@echo "✓ Build complete"

# Run simulation + prepare gen0 data
test: build
	@echo "→ Running simulation and generating per-generator CSV files"
	@./test
	@echo "✓ Simulation complete, per-generator CSV files in sim/ directory"

# Run all analysis scripts
analysis:
	@echo "→ Running core analysis scripts"
	@mkdir -p report_plot
	@cp sim/gen1.csv PLOT_SAVED_BUS.csv  # Use gen1 for analysis by default #we wont use gen0 as its connected to PV Slack bus
	@python3 scripts/vref_step_analysis.py
	@python3 scripts/pm_step_analysis.py
	@python3 scripts/fault_analysis.py
	@python3 scripts/stability_fault_analysis.py
	@python3 scripts/vref_analysis.py
	@python3 scripts/plot_analysis.py
	@python3 scripts/emf_trace.py
	@python3 scripts/trace_signal_path.py
	@python3 scripts/high_precision_plot.py
	@echo "✓ Analysis complete – see report_plot/"

# Clean everything
clean:
	@echo "→ Cleaning all build artifacts and outputs"
	@rm -f test *.o *.png *.csv
	@rm -rf report_plot sim
	@echo "✓ Clean complete" 