.SUFFIXES:
.SUFFIXES: .f90 .o

FC = gfortran
FCFLAGS = -g
EXECUTABLE = ShearDSMC2

MAINDIR = /home/mgbolase/desktop/DSMC_Shear3/model
MODDIR = /home/mgbolase/desktop/DSMC_Shear3/model/mod
OUTDIR = /home/mgbolase/desktop/DSMC_Shear3/model/measure

MODOBJS = \
constants_mod.o \
DSMC_mod.o \
geometry_mod.o \
output_mod.o \
particle_mod.o \
run_param_mod.o \
gmm_cond_mod.o \
collision_tables_mod.o

OBJS = \
collide.o \
index.o \
initialize.o \
integrate_eom.o \
read_input.o

OUTOBJS = \
measure_init.o \
measure_dem.o \
measure_final.o

ALLOBJS = $(OBJS) $(MODOBJS) $(OUTOBJS)

$(EXECUTABLE) : $(MAINDIR)/main.f90 $(ALLOBJS)
	$(FC) $(ALLOBJS) -o $@ $<

# Modules
###########################################################
particle_mod.o : $(MODDIR)/particle_mod.f90
	$(FC) $(FCFLAGS) -c $<
	
run_param_mod.o : $(MODDIR)/run_param_mod.f90
	$(FC) $(FCFLAGS) -c $<
	
DSMC_mod.o : $(MODDIR)/DSMC_mod.f90
	$(FC) $(FCFLAGS) -c $<
	
geometry_mod.o : $(MODDIR)/geometry_mod.f90
	$(FC) $(FCFLAGS) -c $<

constants_mod.o : $(MODDIR)/constants_mod.f90
	$(FC) $(FCFLAGS) -c $<

gmm_cond_mod.o : $(MODDIR)/gmm_cond_mod.f90
	$(FC) $(FCFLAGS) -c $<

collision_tables_mod.o : $(MODDIR)/collision_tables_mod.f90
	$(FC) $(FCFLAGS) -c $<

# Other
###########################################################
output_mod.o : $(MODDIR)/output_mod.f90 $(MODOBJS)
	$(FC) $(FCFLAGS) -c $<

initialize.o : $(MAINDIR)/initialize.f90 $(MODOBJS)
	$(FC) $(FCFLAGS) -c $<

read_input.o : $(MAINDIR)/read_input.f90 $(MODOBJS)
	$(FC) $(FCFLAGS) -c $<
	

# DES
###########################################################
integrate_eom.o : $(MAINDIR)/integrate_eom.f90 $(MODOBJS)
	$(FC) $(FCFLAGS) -c $<

# DSMC
###########################################################
collide.o : $(MAINDIR)/collide.f90 $(MODOBJS)
	$(FC) $(FCFLAGS) -c $<
	
index.o : $(MAINDIR)/index.f90 $(MODOBJS)
	$(FC) $(FCFLAGS) -c $<
	
# Outputs
###########################################################
measure_init.o : $(OUTDIR)/measure_init.f90 $(MODOBJS)
	$(FC) $(FCFLAGS) -c $<
	
measure_dem.o : $(OUTDIR)/measure_dem.f90 $(MODOBJS)
	$(FC) $(FCFLAGS) -c $<
	
measure_final.o : $(OUTDIR)/measure_final.f90 $(MODOBJS)
	$(FC) $(FCFLAGS) -c $<
	
.PHONY : clean

clean :
	rm *.mod *.o