using LinearAlgebra: eigen, Symmetric

struct PODParams
    species::Vector{Symbol}
    rin::Float64
    rcut::Float64
    Pα::Int64
    Pγ::Int64
    Nr2::Int64
    ns::Int64
    # inner constructor will be needed to enforce Nr3 <=Nr2, etc.
    function PODParams(species::Vector{Symbol},
                      rin::Float64,
                      rcut::Float64,
                      Pα::Int64,
                      Pγ::Int64,
                      Nr2::Int64)
        # This is fixed
        Pβ = 3
        # number of snapshots
        ns = Pα*Pβ + Pγ 

        new(species,rin,rcut,Pα,Pγ,Nr2,ns)
    end
end       

struct PODBasis
    params::PODParams
    Φ::Matrix{Float64}

    # inner constructor will be needed to enforce Nr3 <=Nr2, etc.
    function PODBasis(params)
        Φ = init2body(params)
        return new(params,Φ)
   end
end 

function PODBasis(species, rin, rcut, Pα, Pγ, Nr2)
    params = PODParams(species, rin, rcut, Pα, Pγ, Nr2)
    return PODBasis(params)
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
end

function compute_force_descriptors(sys, basis::PODBasis;
                                   neighbors=nothing)
end

# initialization only, can be done on CPU
function init2body(params)

    Φ = eigenvaluedecomposition!(2000,params)
    return Φ   
end

# initialization only, can be done on CPU
function eigenvaluedecomposition!(N,params)
    ns = params.ns    
    rin = params.rin
    rcut = params.rcut

    # `snapshots` assumes S is a contiguous vector 
    S = Vector{Float64}(undef,N*ns)

    xij = [(rin + 1e-6) + (rcut - rin - 1e-6) * (i / (N - 1)) for i in 0:(N-1)]

    snapshots!(S,xij,N, params)
    
    # convert to a matrix for simplicity
    S = reshape(S,N,params.ns)

    # Compute S^TS (i.e. covariance matrix) 
    A = (1/N).*(S'*S)
    
    #Get eigenvectors
    _λ, Φ = eigen(Symmetric(A,:U), sortby=-)

    # project snapshots onto orthogonal basis
    Q = S*Φ
    
    # compute area to normalize the eignevectors
    # should be the same value at every entry?
    xij_diff  = [xij[i+1]-xij[i] for i in 1:N-1]
    for m in 1:ns
        area = 0.0
        for i in 1:N-1
            area += 0.5*xij_diff[i]*(Q[i,m]*Q[i,m] + Q[i+1,m]*Q[i+1,m])
        end
        Φ[:,m] ./= sqrt(area)
    end
   
    # enforce consistent signs for eigenvectors
    for m in 1:ns
        if Φ[m,m] < 0.0
            Φ[:,m] .*= -1
        end
    end

    return Φ
end

# initialization only, can be done on CPU
function snapshots!(S,xij,N, params)
    rcut = params.rcut 
    rin  = params.rin
    rmax = rcut - rin
    
    # βk values for k=1,2,3
    besselparams = [1e-3, 2.0, 4.0]

    Pα = params.Pα
    Pγ = params.Pγ
    
    # iterate through the xij grid
    for n in eachindex(xij)
        dij = xij[n]

        r = dij - rin;

        #compute the normalized distance
        y = r/rmax;

        # cutoff function 
        y2 = y*y;
        y3 = 1.0 - y2*y;
        y4 = y3*y3 + 1e-6;
        y5 = sqrt(y4);
        y6 = exp(-1.0/y5);
        fcut = y6/exp(-1.0);
        
        # iterating over betas, fixed Pβ=3
        for k in 1:3 
            βk = besselparams[k]
            #this is useless?
            if(abs(βk) <= 1e-6)
                βk = 1e-3
            end

            x = (1.0 - exp(-βk*r/rmax))/(1.0-exp(-βk))

            for i in 1:Pα
                a = i*π
                b = (sqrt(2.0/(rmax))/i)

                nij = n + N*(i-1) + N*Pα*(k-1)

                S[nij] = b*fcut*sin(a*x)/r
            end
        end

        for i in 1:Pγ
            p = Pα*3 + (i-1)
            nij = n + N*p
            a = dij^i

            S[nij] = fcut/a
        end
    end
end    
