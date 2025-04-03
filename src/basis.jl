#relevant counts used throughout the calculations
struct PODCounts
    nelements::Int64
    ns::Int64
    nabf3::Int64
    nabf4::Int64
    K3::Int64
    K4::Int64
    nl2::Int64
    nl3::Int64
    nl4::Int64
    nCoeffPerElement::Int64
end

function initialize_counts(species, besseldegree, inversedegree, 
                           nrbf2, nrbf3, nrbf4, P3,P4)
    # These values are fixed
    Pβ = 3
    npa = (0,1,4,10,20,35,56,84,120,165,220,286,364,455)
    nb  = (1,2,4,7,11,16,23)

    # number of snapshots
    ns = besseldegree*Pβ + inversedegree

    #number of elements
    nelements = Ne = length(species)

    # number of angular basis functions
    nabf3 = P3 + 1
    nabf4 = nb[P4+1]

    #number of monomials
    K3 = npa[P3+2]
    K4 = npa[P4+2]

    # number of descriptors
    nl2 = nrbf2*Ne
    nl3 = nabf3*nrbf3*Ne*(Ne+1)÷2
    nl4 = nabf4*nrbf4*Ne*(Ne+1)*(Ne+2)÷6
    nCoeffPerElement = 1 + nl2 + nl3 + nl4

    return PODCounts(nelements,ns,nabf3,nabf4,K3,K4,
                     nl2,nl3,nl4,nCoeffPerElement)

  
end

#parameters of the potentials
struct PODParams
    species::Vector{Symbol}
    rin::Float64
    rcut::Float64
    besseldegree::Int64 #Pα
    inversedegree::Int64 #Pγ
    nrbf2::Int64
    nrbf3::Int64
    nrbf4::Int64
    P3::Int64
    P4::Int64
    counts::PODCounts

    # inner constructor
    function PODParams(species::Vector{Symbol},
                      rin::Float64,
                      rcut::Float64,
                      besseldegree::Int64,
                      inversedegree::Int64,
                      nrbf2::Int64,
                      nrbf3::Int64,
                      nrbf4::Int64,
                      P3::Int64,
                      P4::Int64)

        if nrbf3 > nrbf2 
            error("number of 3-body rbfs has to be <= number of 2-body rbfs")        
        end

        if nrbf4 > nrbf3
            error("number of 3-body rbfs has to be <= number of 2-body rbfs")        
        end

        if P4 > P3 
            error("4-body angular degree must be <= 3-body angular degree")
        end

        if P3 > 12 
            error("3-body angular degree must be <=12")
        end

        if P4 > 6 
            error("4-body angular degree must be <=6")
        end

        counts = initialize_counts(species, besseldegree, inversedegree, 
                                   nrbf2, nrbf3, nrbf4, P3,P4)
        new(species,rin,rcut,besseldegree,inversedegree,
            nrbf2,nrbf3, nrbf4, P3, P4, counts)
    end
end       

include("./pod_initialization.jl")

mutable struct PODWorkspace
    ### size is fixed by params
    bd::Vector{Float64} 

    ### size is a function of nijmax
    bdd::Vector{Float64} # why does this change with Nij?

    # can de-comment when the rbft*phi matmul operation 
    # is non-allocating and directly outputs flat array
    #rbf::Vector{Float64}
    #rbfx::Vector{Float64}
    #rbfy::Vector{Float64}
    #rbfz::Vector{Float64}

    rbft::Vector{Float64}
    rbfxt::Vector{Float64}
    rbfyt::Vector{Float64}
    rbfzt::Vector{Float64}

end

function initialize_workspace(params::PODParams)
    Mdesc = params.counts.nCoeffPerElement -1     
    bd = Vector{Float64}(undef,Mdesc)

    return PODWorkspace(bd,
                        Vector{Float64}(undef,0), #bdd
                        #Vector{Float64}(undef,0), # rbf
                        #Vector{Float64}(undef,0), # rbfx
                        #Vector{Float64}(undef,0), # rbfy
                        #Vector{Float64}(undef,0), # rbfz
                        Vector{Float64}(undef,0), # rbft
                        Vector{Float64}(undef,0), # rbfxt
                        Vector{Float64}(undef,0), # rbfxy
                        Vector{Float64}(undef,0)) # rbfxz
end


struct PODBasis
    params::PODParams
    Φ::Matrix{Float64}
    workspace::PODWorkspace

    # inner constructor will be needed to enforce Nr3 <=Nr2, etc.
    function PODBasis(params)
        Φ = init2body(params)
        workspace = initialize_workspace(params)
        return new(params,Φ,workspace)
   end
end 

function PODBasis(species,rcut; 
                  rin=1.0, 
                  besseldegree=4, 
                  inversedegree=8,
                  nrbf2=8,
                  nrbf3=6,
                  nrbf4=0,
                  P3=4,
                  P4=0)
    params = PODParams(species,rin,rcut,besseldegree, inversedegree,
                       nrbf2, nrbf3, nrbf4, P3, P4)
    return PODBasis(params)
end


function allocate_workspace!(basis::PODBasis,nijmax)
    workspace = basis.workspace

    params = basis.params
    nrbf2 = params.nrbf2
    ns = params.counts.ns
    Mdesc = params.counts.nCoeffPerElement - 1

    workspace.bdd = Vector{Float64}(undef,3*nijmax*Mdesc)

    # can de-comment when the rbft*phi matmul operation 
    # is non-allocating and directly outputs flat array
    #workspace.rbf  = Vector{Float64}(undef,nrbf2*nijmax)
    #workspace.rbfx = Vector{Float64}(undef,nrbf2*nijmax)
    #workspace.rbfy = Vector{Float64}(undef,nrbf2*nijmax)
    #workspace.rbfz = Vector{Float64}(undef,nrbf2*nijmax)

    workspace.rbft  = Vector{Float64}(undef,ns*nijmax)
    workspace.rbfxt = Vector{Float64}(undef,ns*nijmax)
    workspace.rbfyt = Vector{Float64}(undef,ns*nijmax)
    workspace.rbfzt = Vector{Float64}(undef,ns*nijmax)
end

function peratom_length(basis::PODBasis)
    return basis.params.counts.nCoeffPerElement 
end

function Base.length(basis::PODBasis)
    nelements = basis.params.counts.nelements
    nCoeffPerElement = basis.params.counts.nCoeffPerElement
    return nelements*nCoeffPerElement
end

function peratomenergyforce2!(fij, rij,ti,tj,Nj,basis,coeff)
    
    counts = basis.params.counts 
    nCoeffPerElement = counts.nCoeffPerElement
    nelements = counts.nCoeffPerElement
    
    nl2 = counts.nl2
    ns  = counts.ns 
    nrbf2 = basis.params.nrbf2

    workspace = basis.workspace

    # initialize energies and forces to 0.0
    e = 0.0
    for n in 1:3*Nj
        fij[n] = 0.0
    end

    if Nj == 0 
        # just a one-body contribution
        return coeff[nCoeffPerElement*(ti[0]-1)+1] # one-index adjustment
    end
    
    # TODO Allocations that need to be pre-allocated for the GPU

    # This logic only works for the 2-body, will update for 3-body/4-body
    #d2 = bd = Vector{Float64}(undef,nl2)
    #cb2 = cb = bdd = Vector{Float64}(undef, nl2)
    bd = workspace.bd
    d2 = view(bd,1:nl2)

    cb = workspace.bdd
    cb2 = view(cb,1:nl2) # why does bdd have Nj dependence

    #rbf = Vector{Float64}(undef,Nj*basis.params.nrbf2)
    #rbfx = Vector{Float64}(undef,Nj*basis.params.nrbf2)
    #rbfy = Vector{Float64}(undef,Nj*basis.params.nrbf2)
    #rbfz = Vector{Float64}(undef,Nj*basis.params.nrbf2)

    #$rbft  = Vector{Float64}(undef,Nj*ns)
    #$rbfxt = Vector{Float64}(undef,Nj*ns)
    #$rbfyt = Vector{Float64}(undef,Nj*ns)
    #$rbfzt = Vector{Float64}(undef,Nj*ns)

    rbft  = view(workspace.rbft,  1:Nj*ns)
    rbfxt = view(workspace.rbfxt, 1:Nj*ns)
    rbfyt = view(workspace.rbfyt, 1:Nj*ns)
    rbfzt = view(workspace.rbfzt, 1:Nj*ns)

    
    radialbasis!(rbft,rbfxt,rbfyt,rbfzt,rij,basis.params,Nj)
    
    # Apologies, apologies... this will need to be updated for GPU
    rbf =  reshape( reshape(rbft, Nj, ns)*basis.Φ[:, 1:nrbf2], :)
    rbfx = reshape( reshape(rbfxt, Nj, ns)*basis.Φ[:, 1:nrbf2],:)
    rbfy = reshape( reshape(rbfyt, Nj, ns)*basis.Φ[:, 1:nrbf2],:)
    rbfz = reshape( reshape(rbfzt, Nj, ns)*basis.Φ[:, 1:nrbf2],:)

    if nl2> 0 && Nj > 2 
        twobodydesc!(d2, rbf, tj, Nj, basis.params) # extra argument
    end

    e +=  peratombase_coefficients(cb,bd,ti,basis.params,coeff) # two extra arguments

    if nl2> 0 && Nj > 2 
        twobody_forces!(fij,cb2, rbfx, rbfy, rbfz, tj, Nj, nrbf2) # extra arguments
    end

    return e
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

function twobodydesc!(d2,rbf,tj,N, params)
    nelements = params.counts.nelements
    nrbf2 = params.nrbf2
    nl2 = params.counts.nl2
    
    # initialize the two body descriptors 
    for m in 1:nl2
        d2[m] = 0.0
    end

    for m in 1:nrbf2
        for n in 1:N
            i2 = n + N*(m-1) # adjustment for 1-indexing
            d2[m+nrbf2*(tj[n]-1)] += rbf[i2]
        end
   end
end

function peratombase_coefficients(cb,bd,ti,params,coeff)
    nCoeffPerElement = params.counts.nCoeffPerElement
    Mdesc = nCoeffPerElement - 1 # remove one-body term
    nc = nCoeffPerElement*(ti[1]-1)

    ei = coeff[1+nc] 
    for m in 1:Mdesc
        ei += coeff[1 + m + nc]*bd[m]
        cb[m] = coeff[1 + m + nc]
    end
      return ei
end

function twobody_forces!(fij,cb2,rbfx,rbfy,rbfz,tj,Nj,nrbf2)
    totalIterations = nrbf2*Nj
    for idx in 0:totalIterations-1 # keeping this starting at 0 
        n = idx ÷ nrbf2 + 1  # adjust for 1-indexing
        m = mod(idx, nrbf2) + 1 # adjust for 1-indexing

        i2 = n + Nj*(m-1) # n already had the 1-index adjustment
        i1 = 3*(n-1) # don't do 1-indexing adjustmnet here, do it in the fij
        c = cb2[m+nrbf2*(tj[n]-1)] # m already has the 1-index adjustment
        fij[1+i1] += c*rbfx[i2] 
        fij[2+i1] += c*rbfy[i2]
        fij[3+i1] += c*rbfz[i2]
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
