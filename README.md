# 3D MHD simulations of kink waves cutoff through the solar transition region with PLUTO.

This repository contains the setups and analysis code used for the paper [Pelouze et al. 2023, A&A 672, A105](https://doi.org/10.1051/0004-6361/202245049).

This features files to simulate the propagation of kink waves through the transition region, and study their cutoff.

## Directories description

- `setups/`: PLUTO setup files (see their respective README):
  - `open_chromos_relax_2D`: 2D flux tube relaxation
  - `open_chromos_3D`: 3D kink wave propagation
- `conversion/relax_2D_to_3D.py`: converts the output of the 2D relaxation runs to input files for the 3D simulations.
- `mock/`: Python scripts that replicate some C functions from the PLUTO setup.
- `viz/`: Python scripts to analyze simulation output and make plots:
  - `viz/paper_kink_cutoff.py` and `viz/paper_vrw_av.py`: makes the plots shown in the paper
  - `viz/amplitude_altitude_cut_3D.py` computes the wave amplitude (`ampl_alt_data.npz`)
  - `viz/phase_altitude_3D.py` computes the phase `phase_amplitude_data.npz` files
  - remaining files are dependencies of the above.
