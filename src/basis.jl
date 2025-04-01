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

function peratom_length(basis::PODBasis)
    nrbf2 = basis.params.nrbf2 
    nelements = length(basis.params.species)
    
    # one-body and two-body
    return nelements + nrbf2*nelements
end

function Base.length(basis::PODBasis)
    nelements = length(basis.params.species)
    nCoeffPerElement = peratom_length(basis)
    return nelements*nCoeffPerElement
end

function peratomenergyforce2!(fij, rij,ti,tj,Nj,basis,coeff)

    nCoeffPerElement = peratom_length(basis)
    nelements = length(basis.params.species)
    
    nl2 = basis.params.nrbf2 * nelements 


    for n in 1:3*Nj
        fij[n] = 0.0
    end

    if Nj == 0 
        # just a one-body contribution
        return coeff[nCoeffPerElement*(ti[0]-1)+1] # one-index adjustment
    end
    
    # TODO Allocations that need to be pre-allocated for the GPU

    # This logic only works for the 2-body, will update for 3-body/4-body
    d2 = bd = Vector{Float64}(undef,nl2)
    rbf = Vector{Float64}(undef,Nj*basis.params.nrbf2)
    rbfx = Vector{Float64}(undef,Nj*basis.params.nrbf2)
    rbfy = Vector{Float64}(undef,Nj*basis.params.nrbf2)
    rbfz = Vector{Float64}(undef,Nj*basis.params.nrbf2)

    rbft = Vector{Float64}(undef,Nj*basis.params.ns)
    rbfxt = Vector{Float64}(undef,Nj*basis.params.ns)
    rbfyt = Vector{Float64}(undef,Nj*basis.params.ns)
    rbfzt = Vector{Float64}(undef,Nj*basis.params.ns)
    
    radialbasis!(rbft,rbfxt,rbfyt,rbfzt,rij,basis.params,Nj)
end

function radialbasis!(rbf,rbfx,rbfy,rbfz,rij,params,N)

    besselparams = [1e-3, 2.0, 4.0]
    rin = params.rin
    rmax = params.rcut - rin
    besseldegree = params.besseldegree
    inversedegree = params.inversedegree

    nbesselpars = 3

    for n in 1:N 
        xij1 = rij[1+3*(n-1)]
        xij2 = rij[2+3*(n-1)]
        xij3 = rij[3+3*(n-1)]

        dij = sqrt(xij1*xij1 + xij2*xij2 + xij3*xij3)
        dr1 = xij1/dij
        dr2 = xij2/dij
        dr3 = xij3/dij

        r = dij - rin
        y = r/rmax
        y2 = y*y
        
        # compute cutoff function and derivative
        y3 = 1.0 - y2*y
        y4 = y3*y3 + 1e-6
        y5 = sqrt(y4)
        y6 = exp(-1.0/y5)
        y7 = y4*sqrt(y4)
        fcut = y6/exp(-1.0)
        dfcut = ((3.0/(rmax*exp(-1.0)))*(y2)*y6*(y*y2 - 1.0))/y7
        
        # computing the zeroth order bessel function snapshots and derivatives
        f1 = fcut/r
        f2 = f1/r
        df1 = dfcut/r

        alpha = besselparams[1]
        t1 = (1.0-exp(-alpha))
        t2 = exp(-alpha*r/rmax)
        x0 =  (1.0 - t2)/t1
        dx0 = (alpha/rmax)*t2/t1

        alpha = besselparams[2]
        t1 = (1.0-exp(-alpha))
        t2 = exp(-alpha*r/rmax)
        x1 =  (1.0 - t2)/t1
        dx1 = (alpha/rmax)*t2/t1

        alpha = besselparams[3]
        t1 = (1.0-exp(-alpha))
        t2 = exp(-alpha*r/rmax)
        x2 =  (1.0 - t2)/t1
        dx2 = (alpha/rmax)*t2/t1
        for i in 0:besseldegree-1
            a = (i+1)*π
            b = (sqrt(2.0/(rmax))/(i+1))
            af1 = a*f1

            sinax = sin(a*x0)
            nij = n + N*i # this is OK because n starts at 1, i starts at 0

            rbf[nij] = b*f1*sinax
            drbfdr = b*(df1*sinax - f2*sinax + af1*cos(a*x0)*dx0)
            rbfx[nij] = drbfdr*dr1
            rbfy[nij] = drbfdr*dr2
            rbfz[nij] = drbfdr*dr3

            sinax = sin(a*x1)
            nij = n + N*i + N*besseldegree*1 # this is OK because n starts at 1, i starts at 0

            rbf[nij] = b*f1*sinax
            drbfdr = b*(df1*sinax - f2*sinax + af1*cos(a*x1)*dx1)
            rbfx[nij] = drbfdr*dr1
            rbfy[nij] = drbfdr*dr2
            rbfz[nij] = drbfdr*dr3

            sinax = sin(a*x2)
            nij = n + N*i + N*besseldegree*2 # this is OK because n starts at 1, i starts at 0
            rbf[nij] = b*f1*sinax
            drbfdr = b*(df1*sinax - f2*sinax + af1*cos(a*x2)*dx2)
            rbfx[nij] = drbfdr*dr1
            rbfy[nij] = drbfdr*dr2
            rbfz[nij] = drbfdr*dr3
        end

        # computing inverse snapshots and derivatives 
        f1 = fcut/dij
        for i in 0:inversedegree-1
            p = besseldegree*nbesselpars + i
            nij = n + N*p # OK because n starts at 1, i starts at 0
            a = dij^(i+1)

            rbf[nij] = fcut/a

            drbfdr = (dfcut - (i+1.0)*f1)/a
            rbfx[nij] = drbfdr*dr1
            rbfy[nij] = drbfdr*dr2
            rbfz[nij] = drbfdr*dr3        
        end
    end
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
