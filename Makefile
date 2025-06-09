.PHONY: all clean test analysis plots_gen

# Default: build the executable
all: build

# Build directory for object files
OBJ_DIR := obj

# Source files
SOURCES := src/core/main.c src/io/read_fn.c src/math/invert.c src/core/simulator.c \
           src/io/user_input.c src/core/fault_config.c src/network/y_bus_utils.c \
           src/generators/gen_init.c src/network/fault_matrix.c src/generators/power_calc.c \
           src/network/network_solver.c src/generators/current_calc.c

# Object files in build directory
OBJECTS := $(OBJ_DIR)/main.o $(OBJ_DIR)/read_fn.o $(OBJ_DIR)/invert.o $(OBJ_DIR)/simulator.o \
           $(OBJ_DIR)/user_input.o $(OBJ_DIR)/fault_config.o $(OBJ_DIR)/y_bus_utils.o \
           $(OBJ_DIR)/gen_init.o $(OBJ_DIR)/fault_matrix.o $(OBJ_DIR)/power_calc.o \
           $(OBJ_DIR)/network_solver.o $(OBJ_DIR)/current_calc.o

# Compile everything
build: executable

executable: $(OBJECTS)
	@echo "→ Linking executable"
	@mkdir -p report_plot sim
	@gcc $(OBJECTS) -lgsl -lgslcblas -lm -o test
	@echo "✓ Build complete"

# Create build directory
$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

# Individual compilation rules - all output to build directory
$(OBJ_DIR)/main.o: src/core/main.c | $(OBJ_DIR)
	@gcc -g -Wall -c src/core/main.c -o $@

$(OBJ_DIR)/read_fn.o: src/io/read_fn.c | $(OBJ_DIR)
	@gcc -g -Wall -c src/io/read_fn.c -o $@

$(OBJ_DIR)/invert.o: src/math/invert.c | $(OBJ_DIR)
	@gcc -g -Wall -c src/math/invert.c -o $@

$(OBJ_DIR)/simulator.o: src/core/simulator.c | $(OBJ_DIR)
	@gcc -g -Wall -c src/core/simulator.c -o $@

$(OBJ_DIR)/user_input.o: src/io/user_input.c | $(OBJ_DIR)
	@gcc -g -Wall -c src/io/user_input.c -o $@

$(OBJ_DIR)/fault_config.o: src/core/fault_config.c | $(OBJ_DIR)
	@gcc -g -Wall -c src/core/fault_config.c -o $@

$(OBJ_DIR)/y_bus_utils.o: src/network/y_bus_utils.c | $(OBJ_DIR)
	@gcc -g -Wall -c src/network/y_bus_utils.c -o $@

$(OBJ_DIR)/gen_init.o: src/generators/gen_init.c | $(OBJ_DIR)
	@gcc -g -Wall -c src/generators/gen_init.c -o $@

$(OBJ_DIR)/fault_matrix.o: src/network/fault_matrix.c | $(OBJ_DIR)
	@gcc -g -Wall -c src/network/fault_matrix.c -o $@

$(OBJ_DIR)/power_calc.o: src/generators/power_calc.c | $(OBJ_DIR)
	@gcc -g -Wall -c src/generators/power_calc.c -o $@

$(OBJ_DIR)/network_solver.o: src/network/network_solver.c | $(OBJ_DIR)
	@gcc -g -Wall -c src/network/network_solver.c -o $@

$(OBJ_DIR)/current_calc.o: src/generators/current_calc.c | $(OBJ_DIR)
	@gcc -g -Wall -c src/generators/current_calc.c -o $@

# Run simulation + prepare gen0 data
test: executable
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
	@rm -f test *.csv
	@rm -rf $(OBJ_DIR) report_plot sim
	@echo "✓ Clean complete" 