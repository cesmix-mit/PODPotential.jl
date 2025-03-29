using AtomsCalculators: potential_energy, 
                        forces, 
                        virial,
                        energy_unit,
                        length_unit

struct PODPotential
    basis::PODBasis
    β::Vector{Float64}
end

energy_unit(calc::PODPotential) = u"eV"
length_unit(calc::PODPotential) = u"Å"

function potential_energy(sys, calc::PODPotential; 
                          neighbors=nothing,
                          kwargs...) 

    gd = compute_global_descriptors(sys,calc;neighbors=neighbors)
    pe = gd ⋅ calc.β
end

function forces(sys, calc::PODPotential;
                neighbors=nothing,
                kwargs...)

    fd = compute_force_descriptors(sys,calc;neighbors=neighbors)  
    return nothing
end

function virial(sys,calc::PODPotential;
                neighbors=nothing
                kwargs...)
    error("virial isn't supported for PODPotential yet")
end
