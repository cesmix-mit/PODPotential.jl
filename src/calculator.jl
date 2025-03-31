import AtomsCalculators: potential_energy, 
                        forces, 
                        virial,
                        energy_unit,
                        length_unit

struct POD_Potential
    basis::PODBasis
    β::Vector{Float64}
end

energy_unit(calc::POD_Potential) = u"eV"
length_unit(calc::POD_Potential) = u"Å"

function potential_energy(sys, calc::POD_Potential; 
                          neighbors=nothing,
                          kwargs...) 

    gd = compute_global_descriptors(sys,calc;neighbors=neighbors)
    pe = gd ⋅ calc.β
end

function forces(sys, calc::POD_Potential;
                neighbors=nothing,
                kwargs...)

    fd = compute_force_descriptors(sys,calc;neighbors=neighbors)  
    return nothing
end

function virial(sys,calc::POD_Potential;
                neighbors=nothing,
                kwargs...)
    error("virial isn't supported for PODPotential yet")
end
