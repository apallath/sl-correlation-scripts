Format `<variable>[whitespace]<value>`

- `snapshots`: number of snapshots in dump file

- `natoms`: number of atoms
- `nliquidatoms`: number of liquid atoms
- `nsolidatoms`: number of solid atoms
- `nsolidtopatoms`: number of atoms in topmost solid layer

- `x_min`: min x value of box
- `x_max`: max x value of box
- `y_min`: min y value of box
- `y_max`: max y value of box
- `z_min`: min z value of box
- `z_max`: max z value of box

- `positions_file`: trajectory file containing atom positions
- `solid_top_positions_file`: trajectory file containing topmost layer solid positions
- `solid_init_top_positions_file`: file containing snapshot of topmost layer initial solid positions
- `msd_file`: dump file containing mean square displacements

- `equil_dump_file_pattern`: naming pattern of equilibration trajectory dump file
- `prod_dump_file_pattern`: naming pattern of production trajectory dump file
