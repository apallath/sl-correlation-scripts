# Input data
(From production run of LAMMPS MD simulations for solid)

- Mean Square Displacement dump files: `MSD_all_{}.dat`
- Topmost solid layer trajectories - scaled coordinates: `toppositions_init_{}.atom`
- Topmost solid layer initial snapshot (before equilibration) - scaled coordinates: `toppositions_{}.lammpstrj`
- Solid-liquid system trajectories - scaled coordinates: `positions_{}.atom`

# Correlations
- displacement-displacement correlations for topmost solid layer
- z_min - z_min correlations of liquid above topmost solid layer surface 

# Run all
Run the MATLAB file `driver.m`



