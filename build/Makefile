# DSMC_NSP Makefile with portable relative paths
# Build from build/ directory: make
# Clean: make clean
# Remove everything: make cleanall

FC = gfortran
FFLAGS = -O3 -march=native -ffast-math -funroll-loops

# Directory structure (relative to build/)
SRCDIR = ../src
MAINDIR = $(SRCDIR)/main
MODDIR = $(SRCDIR)/modules
MEASDIR = $(SRCDIR)/measure
BINDIR = ../bin
BUILDDIR = .

# Module output directory
MODFLAG = -J$(BUILDDIR)

# Target executable
TARGET = $(BINDIR)/DSMC

# Module sources (must be compiled first, in dependency order)
MOD_SOURCES = \
	$(MODDIR)/constants_mod.f90 \
	$(MODDIR)/geometry_mod.f90 \
	$(MODDIR)/particle_mod.f90 \
	$(MODDIR)/run_param_mod.f90 \
	$(MODDIR)/DSMC_mod.f90 \
	$(MODDIR)/gmm_cond_mod.f90 \
	$(MODDIR)/collision_tables_mod.f90 \
	$(MODDIR)/output_mod.f90

# Main sources
MAIN_SOURCES = \
	$(MAINDIR)/read_input.f90 \
	$(MAINDIR)/initialize.f90 \
	$(MAINDIR)/index.f90 \
	$(MAINDIR)/collide.f90 \
	$(MAINDIR)/integrate_eom.f90

# Measurement sources
MEAS_SOURCES = \
	$(MEASDIR)/measure_init.f90 \
	$(MEASDIR)/measure_dem.f90 \
	$(MEASDIR)/measure_final.f90

# Main program
PROG_SOURCE = $(MAINDIR)/main.f90

# Object files
MOD_OBJS = $(notdir $(MOD_SOURCES:.f90=.o))
MAIN_OBJS = $(notdir $(MAIN_SOURCES:.f90=.o))
MEAS_OBJS = $(notdir $(MEAS_SOURCES:.f90=.o))
PROG_OBJ = main.o

ALL_OBJS = $(MOD_OBJS) $(MAIN_OBJS) $(MEAS_OBJS) $(PROG_OBJ)

# Default target
all: $(TARGET)

# Link executable
$(TARGET): $(ALL_OBJS) | $(BINDIR)
	$(FC) $(FFLAGS) -o $@ $(ALL_OBJS)
	@echo "Build complete: $@"

# Create bin directory if needed
$(BINDIR):
	mkdir -p $(BINDIR)

# Compile modules (order matters!)
constants_mod.o: $(MODDIR)/constants_mod.f90
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

geometry_mod.o: $(MODDIR)/geometry_mod.f90 constants_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

particle_mod.o: $(MODDIR)/particle_mod.f90 constants_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

run_param_mod.o: $(MODDIR)/run_param_mod.f90
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

DSMC_mod.o: $(MODDIR)/DSMC_mod.f90 constants_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

gmm_cond_mod.o: $(MODDIR)/gmm_cond_mod.f90
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

collision_tables_mod.o: $(MODDIR)/collision_tables_mod.f90 constants_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

output_mod.o: $(MODDIR)/output_mod.f90 constants_mod.o particle_mod.o geometry_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

# Main sources
read_input.o: $(MAINDIR)/read_input.f90 geometry_mod.o particle_mod.o DSMC_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

initialize.o: $(MAINDIR)/initialize.f90 particle_mod.o constants_mod.o DSMC_mod.o geometry_mod.o collision_tables_mod.o gmm_cond_mod.o run_param_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

index.o: $(MAINDIR)/index.f90 particle_mod.o DSMC_mod.o geometry_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

collide.o: $(MAINDIR)/collide.f90 particle_mod.o DSMC_mod.o constants_mod.o geometry_mod.o run_param_mod.o collision_tables_mod.o gmm_cond_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

integrate_eom.o: $(MAINDIR)/integrate_eom.f90 particle_mod.o geometry_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

# Measurement sources
measure_init.o: $(MEASDIR)/measure_init.f90 output_mod.o particle_mod.o geometry_mod.o DSMC_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

measure_dem.o: $(MEASDIR)/measure_dem.f90 output_mod.o particle_mod.o geometry_mod.o DSMC_mod.o run_param_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

measure_final.o: $(MEASDIR)/measure_final.f90 output_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

# Main program
main.o: $(MAINDIR)/main.f90 run_param_mod.o
	$(FC) $(FFLAGS) $(MODFLAG) -c $<

# Clean targets
clean:
	rm -f *.o *.mod

cleanall: clean
	rm -f $(TARGET)

.PHONY: all clean cleanall
