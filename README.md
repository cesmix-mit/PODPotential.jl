A fairly direct port of the ML-POD potential originally written in C++ in LAMMPS.
I try to follow the logic of that code fairly directly, even at the expense of doing things in non-optimal or non-Julian ways.

Restricted to 3-body currently (likely will only consider 4-body max, ignoring the quadratic descriptors)
Not considering environment-dependent descriptors or Gaussian snapshots.

For now, not using the InteratomicPotentials.jl interface, but will eventually (after some updates to IP.jl to handle neighbors)
Satisfies the most basic version of the AtomsCalculators.jl interface, to be compatible with Molly.jl
