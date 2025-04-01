struct PODParams
    species::Vector{Symbol}
    rin::Float64
    rcut::Float64
    besseldegree::Int64 #Pα
    inversedegree::Int64 #Pγ
    nrbf2::Int64
    ns::Int64
    # inner constructor will be needed to enforce Nr3 <=Nr2, etc.
    function PODParams(species::Vector{Symbol},
                      rin::Float64,
                      rcut::Float64,
                      besseldegree::Int64,
                      inversedegree::Int64,
                      nrbf2::Int64)
        # This is fixed
        Pβ = 3
        # number of snapshots
        ns = besseldegree*Pβ + inversedegree

        new(species,rin,rcut,besseldegree,inversedegree,nrbf2,ns)
    end
end       

include("./pod_initialization.jl")

struct PODBasis
    params::PODParams
    Φ::Matrix{Float64}

    # inner constructor will be needed to enforce Nr3 <=Nr2, etc.
    function PODBasis(params)
        Φ = init2body(params)
        return new(params,Φ)
   end
end 

function PODBasis(species, rin, rcut, besseldegree, inversedegree, Nr2)
    params = PODParams(species, rin, rcut, besseldegree, inversedegree, Nr2)
    return PODBasis(params)
end


struct PODWorkspace
end

function Base.length(basis::PODBasis)
    nrbf2 = basis.params.nrbf2
    nelements = length(basis.params.species)

    # only valid for 2-body atm
    return nrbf2*nelements^2
end

### These functions below will be preferred later. For now, just stubs
function compute_local_descriptors(sys, basis::PODBasis;
                                   neighbors=nothing)
end

function compute_global_descriptors(sys, basis::PODBasis;
                                   neighbors=nothing)
end

function compute_force_descriptors(sys, basis::PODBasis;
                                   neighbors=nothing)
end
