
struct PODBasis
end 

struct PODWorkspace
end

function Base.length(basis::PODBasis)
end

function compute_local_descriptors(sys, basis::PODBasis;
                                   neighbors=nothing)
end

function compute_global_descriptors(sys, basis::PODBasis;
                                   neighbors=nothing)

    ld = compute_local_descriptors(sys,basis;neighbors=neighbors)
    gd = sum(ld) # This won't be right when ld is a matrix 
    gd 
end

function compute_force_descriptors(sys, basis::PODBasis;
                                   neighbors=nothing)
end
