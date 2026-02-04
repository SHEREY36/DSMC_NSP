# DSMC-NSP: Direct Simulation Monte Carlo for Non-Spherical Particles

A Fortran 90 simulation framework for granular shear flows of non-spherical particles, coupling DSMC collision dynamics with DEM particle integration using pre-trained Gaussian Mixture Models (GMMs) for collision outcome sampling.

## Directory Structure

```
DSMC_NSP/
├── src/                      # All source code
│   ├── main/                 # Main simulation routines
│   │   ├── main.f90          # Program entry point
│   │   ├── initialize.f90    # System initialization
│   │   ├── read_input.f90    # Input file parsing
│   │   ├── collide.f90       # DSMC collision detection & processing
│   │   ├── index.f90         # Particle-to-cell assignment
│   │   └── integrate_eom.f90 # DEM integration (velocity Verlet)
│   ├── modules/              # Fortran modules
│   │   ├── constants_mod.f90
│   │   ├── geometry_mod.f90
│   │   ├── particle_mod.f90
│   │   ├── DSMC_mod.f90
│   │   ├── run_param_mod.f90
│   │   ├── gmm_cond_mod.f90         # GMM conditional sampling
│   │   ├── collision_tables_mod.f90 # p(chi), P1hit, delta_eps_max tables
│   │   └── output_mod.f90
│   ├── measure/              # Measurement routines
│   │   ├── measure_init.f90
│   │   ├── measure_dem.f90
│   │   └── measure_final.f90
│   └── particle_gen/         # Particle initial condition generator
│       ├── DESParticleGen_MP.f90
│       ├── particles.f90
│       └── test_inputs.f90
├── collision_models/         # Pre-trained GMM binary files
│   ├── gmm_cond_AR11.bin     # AR = 1.1
│   ├── gmm_cond_AR125.bin    # AR = 1.25
│   ├── gmm_cond_AR15.bin     # AR = 1.5
│   ├── gmm_cond_AR20.bin     # AR = 2.0
│   ├── gmm_cond_AR25.bin     # AR = 2.5
│   └── gmm_cond_AR30.bin     # AR = 3.0
├── build/                    # Build directory
│   └── Makefile              # Portable Makefile
├── bin/                      # Compiled executable (created by make)
│   └── DSMC
├── runs/                     # Run directories
│   └── template/             # Template input file
│       └── system_input.dat
├── tools/                    # Helper scripts
│   └── setup_run.sh          # Creates new run directories
├── DSMC_Shear3/              # Legacy structure (preserved for reference)
└── README.md
```

## Quick Start

### 1. Build the Simulation

```bash
cd build
make
```

This creates the executable at `bin/DSMC`.

### 2. Set Up a New Run

```bash
tools/setup_run.sh my_case_name
```

This creates `runs/my_case_name/` with:
- Template `system_input.dat` (edit this for your parameters)
- Symlink to executable (`DSMC`)
- Symlink to collision models

### 3. Configure Your Simulation

Edit `runs/my_case_name/system_input.dat`:

```fortran
# Domain dimensions (m)
0.01 0.01 0.01

# DSMC grid (nx, ny, nz)
10 10 10

# Measurement grid (nx, ny, nz)
10 10 10

# Number of simulated particles
1000

# Scaling factor F_N
1.0

# Particle diameter (m), aspect ratio
0.001 2.0

# Density (kg/m³)
2500.0

# Restitution coefficient
0.9

# Timestep (s), end time (s)
1.0e-5 0.1

# Periodic boundaries (x, y, z): 1=periodic, 0=walls
1 1 1

# Wall velocities (3 lines for x/y/z planes)
0.0 0.0 0.0  0.0 0.0 0.0
0.0 0.0 0.0  0.0 0.0 0.0
0.0 0.0 0.0  0.0 0.0 0.0

# Random walk flag
0
```

**Supported Aspect Ratios:** 1.1, 1.25, 1.5, 2.0, 2.5, 3.0

The correct GMM file is **automatically loaded** based on your AR value.

### 4. Generate Particle Initial Conditions

You need a `particle_input.dat` file with particle positions, velocities, and angular velocities. Use the particle generator or create your own.

### 5. Run the Simulation

```bash
cd runs/my_case_name
./DSMC
```

## Output Files

Each run produces:
- `nu.txt` - Collision frequency
- `tg.txt` - Granular temperature
- `Pijk.txt`, `Pijc.txt` - Pressure tensors (kinetic, collisional)
- `vprof.txt` - Velocity profiles
- `nprof.txt` - Number density profiles
- `tgprof.txt` - Temperature profiles
- `VACF.txt` - Velocity autocorrelation function
- `KIso.txt`, `ICS.txt`, `VoV0.txt` - Additional diagnostics

## Advanced Configuration

### Custom GMM Directory

Set the `DSMC_GMM_DIR` environment variable to use GMM files from a different location:

```bash
export DSMC_GMM_DIR=/path/to/my/gmm/files
./DSMC
```

The code will search for `gmm_cond_AR{label}.bin` in this directory.

### Makefile Targets

```bash
make         # Build executable
make clean   # Remove object files and modules
make cleanall # Remove everything including executable
```

## Migration from Old Structure

If you have existing runs in `DSMC_Shear3/`, you can migrate them:

1. Copy the run directory to `runs/legacy_[name]/`
2. Update the executable symlink:
   ```bash
   cd runs/legacy_[name]
   rm ShearDSMC2  # old executable
   ln -sf ../../bin/DSMC ./DSMC
   ln -sf ../../collision_models ./collision_models
   ```
3. Run with the new executable: `./DSMC`

## Implementation Details

### Collision Model

- **GMM Sampling:** Pre-trained Gaussian Mixture Models map `(r, e_tr, e_r1) → (eps_tr_f, eps_r1_f)` for collision energy redistribution
- **Dense-Regime Correction:** Enskog `chi_ma(solid_fraction)` factor
- **Rejection Sampling:** Uses `p(chi)` probability distribution and `P1hit`, `delta_eps_max` lookup tables
- **Transform:** Logit/sigmoid for bounded [0,1] sampling space

### Boundary Conditions

- **Periodic:** Lees-Edwards with velocity offsets for simple shear
- **Reflective Walls:** Specular reflection with configurable wall velocities

### Particle Model

- **Shape:** Spherocylinders (capped cylinders)
- **Dynamics:** Translational + rotational motion
- **Integration:** Velocity Verlet

### No External Dependencies

Uses only Fortran intrinsics (`random_number`, `dot_product`, etc.). No external libraries required.

## Validation

A standalone validation program tests GMM and collision table functionality:

```bash
cd DSMC_Shear3/model/mod/
gfortran -g -c gmm_cond_mod.f90
gfortran -g -c collision_tables_mod.f90
gfortran -g gmm_cond_mod.o collision_tables_mod.o validate_collision_models.f90 -o validate
./validate
```

Produces: `chi_samples.txt`, `pchi_curve.txt`, `p1hit_table.txt`, `delta_eps_max_table.txt`

## Key Parameters

- **F_N:** Simulation-to-real scaling (real particles = PIP × F_N)
- **AR:** Aspect ratio (length/diameter of spherocylinder)
- **alpha_pp:** Particle-particle restitution coefficient [-1, 1]
- **DSMC grid:** Collision detection resolution (finer = more accurate, slower)
- **Measurement grid:** Output statistics resolution (independent of DSMC grid)

## Troubleshooting

### "GMM file not found"
- Check that AR in `system_input.dat` matches a supported value (1.1, 1.25, 1.5, 2.0, 2.5, 3.0)
- Verify `collision_models/` symlink exists in your run directory
- Check that corresponding `.bin` file exists in `collision_models/`

### "ERROR: cannot open GMM binary file"
- File may be corrupted; re-copy from `DSMC_Shear3/model/mod/`
- Check file permissions

### Build errors
- Ensure you're building from `build/` directory
- Module dependency errors: run `make clean` then `make`

## Citation

If you use this code, please cite the relevant papers (add your citations here).

## License

(Add your license information)

## Contact

(Add contact information)
