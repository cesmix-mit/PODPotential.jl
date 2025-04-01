module PODPotential

using AtomsBase
using Unitful

export LAMMPS_State, lammps_compute
export PODBasis 

include("basis.jl")
include("calculator.jl")

end # module PODPotential
