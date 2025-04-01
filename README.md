Directly replicating LAMMPS logic of the ML-POD potential.

Goal is 3-body, at most 4-body. Ignoring the quadratic descriptors.

Not considering environment-dependent descriptors or Gaussian snapshots.

For now, not using the InteratomicPotentials.jl interface, but will eventually (after some updates to IP.jl to handle neighbors).

Will need to satisfy the most basic version of the AtomsCalculators.jl interface to be compatible with Molly.jl. 

However, intial focus is to take arrays dumped from lammps and directly pass them to routines ported here, to facilitate GPU performance comparison.

TODO:
- Pre-allocate
