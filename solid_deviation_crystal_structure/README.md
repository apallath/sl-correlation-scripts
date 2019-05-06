# Input data
(From production run of LAMMPS MD simulations for solid)

- Mean Square Displacement dump files: `MSD_all_{}.dat`
- Equilibration trajectories: `eq_{}.lammpstrj`
- Production trajectories: `prod_{}.lammpstrj`

# Plot mean square displacement v/s k
- In `plot_msd.py`, enter comma-separated values of `k` in the variable `klist`
- Run `plot_msd.py` using Python 2.7 to generate PNG and PDF plots

# Plot deviation from crystal structure v/s k
- In `crystal_deviation.m`, enter comma-separated values of `k` in the variable `kvals`
- In `crystal_deviation.input`, enter values of variables with comments
- Run `crystal_deviation.m` in MATLAB
