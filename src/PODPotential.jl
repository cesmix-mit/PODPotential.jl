module PODPotential

using AtomsBase
using Unitful

export PODBasis 

include("basis.jl")
include("calculator.jl")

end # module PODPotential
