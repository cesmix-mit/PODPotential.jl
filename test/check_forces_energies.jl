using PODPotential

using JLD2 
using AtomsIO

ref_sys = load_system("./files/reference_monoclinic_hfo2.xyz")

ref_pe = ref_sys.system_data.energy
ref_forces = reduce(hcat, ref_sys.atom_data.forces)'

podbasis = PODBasis([:Hf,:O], 1.0, 5.5, 4,8,8)
podpot = POD_Potential(podbasis,"./files/simple_2body_HfO2_coefficients.pod")
lammps_state = LAMMPS_State("./files/monoclinic_hfo2_lammps_state.jld2")

pe, f = lammps_compute(lammps_state,podpot)
