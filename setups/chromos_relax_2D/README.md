# Overview

In this setup, we derive the pressure and density from the temperature profile
by numerically solving the hydrostatic equilibrium equation.

However, this step cannot easily be parallelized by PLUTO. The differential
equation is integrated from the loop footpoint to the apex (using chromospheric
density as a boundary condition). If PLUTO decomposes the domain along the
loop, the hydrostatic equilibrium will be integrated with a wrong boundary
condition for blocks others than the one containing the footpoint.

To avoid messing with the parallel decomposition, we run two instances of
PLUTO: 

1. Solve the hydrostatic equilibrium by running PLUTO on a single core and save
   the result. We only run the initialization (using `./pluto -maxsteps 0`).
2. Load the output of the previous run and solve the time dependant problem in
   parallel.

To do the above, first run `./setup_init.py`, then `mpicc ./pluto`. 
(No need to run PLUTO’s `setup.py` or make.)

# `setup_init.py`

The `setup_init.py` script generates the required files, compiles PLUTO, and
runs the (non-parallel) hydrostatic equilibrium resolution. It reads
configuration from `setup_config.yml` and performs the following steps:

1. Generate and compile the sequential version of PLUTO (which solves the
   hydrostatic equilibrium) using the `init_hs_equil` section of
   `setup_config.yml`:

    - Generate `pluto.ini` from `pluto.ini.j2` by interpolating the variables
      defined in the `pluto_ini` subsection of `setup_config.yml`.
    - Update `definitions.h` with the variables defined in the `definitions`
      subsection of `setup_config.yml`.
    - (Create output and log directories.)
    - Run PLUTO’s `setup.py` to generate the appropriate makefile and update
      `definitions.h`. (Requires the `$PLUTO_DIR` environement variable.)
    - Run `make` to compile PLUTO.
    - (Run `make clean`.)

2. Run the sequential instance of PLUTO (`./pluto -maxsteps 0`).

3. Generate and compile the parallel version of PLUTO (which solves the time
   dependant equations) using the `solve` section of `setup_config.yml`.
   This is done using the same steps as described in 1.

**Note:** `setup_init.py` **does not** run the parallel version of PLUTO, you
need to later run, eg., `mpirun -np 8 ./pluto`.
