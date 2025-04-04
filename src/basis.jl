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
    ### will never change 
    pn3::Vector{Int64}
    pq3::Vector{Int64}
    pc3::Vector{Int64}
    
    ### size is fixed by params
    bd::Vector{Float64} 
    sumU::Vector{Float64}
    tm::Vector{Float64}
    forcecoeff::Vector{Float64}
    elemindex::Vector{Int64}

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

    U::Vector{Float64}
    Ux::Vector{Float64}
    Uy::Vector{Float64}
    Uz::Vector{Float64}

    abf::Vector{Float64}
    abfx::Vector{Float64}
    abfy::Vector{Float64}
    abfz::Vector{Float64}
end

function initialize_workspace(params::PODParams)
    Mdesc = params.counts.nCoeffPerElement -1     
    K3 = params.counts.K3
    nrbf3 = params.nrbf3
    nelements = Ne = params.counts.nelements

    pn3, pq3, pc3 = init3body(params)

    bd = Vector{Float64}(undef,Mdesc)
    sumU = Vector{Float64}(undef,K3*nrbf3*nelements)
    tm   = Vector{Float64}(undef,4*K3)
    forcecoeff = Vector{Float64}(undef,nelements*K3*nrbf3)

    # keeping elemindex 0-indexed because of how it's used downstream
    elemindex = Vector{Int64}(undef,nelements*nelements)
    k=0
    for i1 in 1:nelements
        for i2 in i1:nelements
            elemindex[i2 + Ne*(i1-1)] = k
            elemindex[i1 + Ne*(i2-1)] = k
            k += 1
        end
    end

    return PODWorkspace(pn3,
                        pq3,
                        pc3,
                        bd,
                        sumU,
                        tm,
                        forcecoeff,
                        elemindex,
                        Vector{Float64}(undef,0), #bdd
                        #Vector{Float64}(undef,0), # rbf
                        #Vector{Float64}(undef,0), # rbfx
                        #Vector{Float64}(undef,0), # rbfy
                        #Vector{Float64}(undef,0), # rbfz
                        Vector{Float64}(undef,0), # rbft
                        Vector{Float64}(undef,0), # rbfxt
                        Vector{Float64}(undef,0), # rbfxy
                        Vector{Float64}(undef,0), # rbfxz
                        Vector{Float64}(undef,0), # U 
                        Vector{Float64}(undef,0), # Ux 
                        Vector{Float64}(undef,0), # Uy 
                        Vector{Float64}(undef,0), # Uz 
                        Vector{Float64}(undef,0), # abf 
                        Vector{Float64}(undef,0), # abfx 
                        Vector{Float64}(undef,0), # abfy 
                        Vector{Float64}(undef,0)) # abfz 
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
    nrbf3 = params.nrbf3 
    K3 = params.counts.K3
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
    
    # pretty sure nl3 > nl4 always...
    if params.counts.nl3 > 0 || params.counts.nl4 > 0
        workspace.U    = Vector{Float64}(undef,K3*nrbf3*nijmax)
        workspace.Ux   = Vector{Float64}(undef,K3*nrbf3*nijmax)
        workspace.Uy   = Vector{Float64}(undef,K3*nrbf3*nijmax)
        workspace.Uz   = Vector{Float64}(undef,K3*nrbf3*nijmax)

        workspace.abf  = Vector{Float64}(undef,K3*nijmax)
        workspace.abfx = Vector{Float64}(undef,K3*nijmax)
        workspace.abfy = Vector{Float64}(undef,K3*nijmax)
        workspace.abfz = Vector{Float64}(undef,K3*nijmax)
    end
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
    nelements = counts.nelements
    
    nl2 = counts.nl2
    nl3 = counts.nl3
    ns  = counts.ns 
    nrbf2 = basis.params.nrbf2
    nrbf3 = basis.params.nrbf3 
    K3 = counts.K3

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
    
    bd = workspace.bd
    d2 = view(bd,1:nl2)

    cb = workspace.bdd
    cb2 = view(cb,1:nl2) # why does bdd have Nj dependence

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

    if nl2 > 0 && Nj > 0 
        twobodydesc!(d2, rbf, tj, Nj, basis.params) # extra argument
    end

    if nl3 > 0  && Nj > 1
        U  = view(workspace.U, 1:Nj*K3*nrbf3)
        Ux = view(workspace.Ux, 1:Nj*K3*nrbf3)
        Uy = view(workspace.Uy, 1:Nj*K3*nrbf3)
        Uz = view(workspace.Uz, 1:Nj*K3*nrbf3)
        sumU = workspace.sumU 

        abf = view(workspace.abf, 1:Nj*K3)
        abfx = view(workspace.abfx, 1:Nj*K3)
        abfy = view(workspace.abfy, 1:Nj*K3)
        abfz = view(workspace.abfz, 1:Nj*K3)
        tm = workspace.tm

        d3 = view(bd,nl2+1:nl2+nl3)
        cb3 = view(cb,nl2+1:nl2+nl3)

        angularbasis!(abf, abfx, abfy, abfz, rij, tm, workspace.pq3, Nj, K3)

        radialangularbasis!(sumU,U,Ux,Uy,Uz,rbf,rbfx,rbfy,rbfz,
                            abf,abfx,abfy,abfz,tj,Nj,K3,nrbf3,nelements)

        threebodydesc!(d3,sumU, workspace.pn3, workspace.pc3, basis.params)
    end

    e +=  peratombase_coefficients(cb,bd,ti,basis.params,coeff) # two extra arguments

    if nl2 > 0 && Nj > 2 
        twobody_forces!(fij,cb2, rbfx, rbfy, rbfz, tj, Nj, nrbf2) # extra arguments
    end

    if nl3 > 0 && Nj > 1
        for i in eachindex(workspace.forcecoeff)
            workspace.forcecoeff[i] = 0.0
        end
        threebody_forcecoeff!(workspace.forcecoeff,cb3,sumU,workspace.pn3, workspace.pc3, workspace.elemindex, basis.params)
        allbody_forces!(fij,workspace.forcecoeff,Ux,Uy,Uz,tj,Nj,basis.params)
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

function angularbasis!(abf, abfx, abfy, abfz, rij, tm, pq, N, K)

    tm0 = view(tm,1:K)
    tmu = view(tm,K+1:2K)
    tmv = view(tm,2K+1:3K)
    tmw = view(tm,3K+1:4K)
    
    # Initialize first angular basis function and its derivatives
    tm0[1] = 1.0
    tmu[1] = 0.0
    tmv[1] = 0.0
    tmw[1] = 0.0

    # loop over all neighboring atoms
    for j in 1:N
        #Calculate relative positions of neighboring atoms and atom i
        x = rij[1+3*(j-1)]
        y = rij[2+3*(j-1)]
        z = rij[3+3*(j-1)]

        # Calculate various terms for derivatives
        xx = x*x
        yy = y*y
        zz = z*z
        xy = x*y
        xz = x*z
        yz = y*z

        # Calculate distance between neighboring atoms and unit vectors
        dij = sqrt(xx + yy + zz)
        u = x/dij
        v = y/dij
        w = z/dij

        # Calculate derivatives of unit vectors
        dij3 = dij*dij*dij
        dudx = (yy+zz)/dij3
        dudy = -xy/dij3
        dudz = -xz/dij3

        dvdx = -xy/dij3
        dvdy = (xx+zz)/dij3
        dvdz = -yz/dij3

        dwdx = -xz/dij3
        dwdy = -yz/dij3
        dwdz = (xx+yy)/dij3

        #Initialize first angular basis function and its derivatives
        abf[j] = tm0[1]
        abfx[j] = 0.0
        abfy[j] = 0.0
        abfz[j] = 0.0

        for n in 2:K # was n=1, n<K
            m = pq[n] # this is an index. Was pq[n] -1
            d = pq[n+K] # this is used as a switch flag

            if d==1
                tm[n] = tm[m]*u
                tmu[n] = tmu[m]*u + tm[m]
                tmv[n] = tmv[m]*u
                tmw[n] = tmw[m]*u
            end

            if d==2
                tm[n] = tm[m]*v
                tmu[n] = tmu[m]*v
                tmv[n] = tmv[m]*v + tm[m]
                tmw[n] = tmw[m]*v
            end
            
            if d==3
                tm[n] = tm[m]*w
                tmu[n] = tmu[m]*w
                tmv[n] = tmv[m]*w
                tmw[n] = tmw[m]*w + tm[m]
            end

            abf[j + N*(n-1)] = tm[n]
            abfx[j + N*(n-1)] = tmu[n]*dudx + tmv[n]*dvdx + tmw[n]*dwdx
            abfy[j + N*(n-1)] = tmu[n]*dudy + tmv[n]*dvdy + tmw[n]*dwdy
            abfz[j + N*(n-1)] = tmu[n]*dudz + tmv[n]*dvdz + tmw[n]*dwdz
        end
    end
end

function radialangularbasis!(sumU,U,Ux,Uy,Uz,rbf,rbfx,rbfy,rbfz,
                             abf,abfx,abfy,abfz,atomtype,N,K,M,Ne)
    for i in eachindex(sumU)
        sumU[i] = 0.0
    end

    #Calculate radial-angular basis functions
    if Ne == 1 # 1-element case (TODO: test)
        for m in 1:M
            for k in 1:K
                sum = 0.0
                for n in 1:N
                    ia = n + N*(k-1)
                    ib = n + N*(m-1)
                    ii = ia + N*K*(m-1)

                    #Calculate c1 and c2
                    c1 = rbf[ib]
                    c2 = abf[ia]

                    #Calculate U, Ux, Uy, Uz
                    U[ii] = c1 * c2
                    Ux[ii] = abfx[ia] * c1 + c2 * rbfx[ib]
                    Uy[ii] = abfy[ia] * c1 + c2 * rbfy[ib]
                    Uz[ii] = abfz[ia] * c1 + c2 * rbfz[ib]

                    #Update sum
                    sum += c1*c2
                end
                #Update sumU
                sumU[k+K*(m-1)] += sum
            end
        end
    else # for multi-element case
        for m in 1:M
            for k in 1:K
                for n in 1:N
                    ia = n + N*(k-1)
                    ib = n + N*(m-1)
                    ii = ia + N*K*(m-1)

                    # Calculate c1 and c2
                    c1 = rbf[ib]
                    c2 = abf[ia]

                    # Calculate U, Ux, Uy, Uz
                    U[ii] = c1*c2
                    Ux[ii] = abfx[ia]*c1 + c2*rbfx[ib]
                    Uy[ii] = abfy[ia]*c1 + c2*rbfy[ib]
                    Uz[ii] = abfz[ia]*c1 + c2*rbfz[ib]

                    # Update sumU with atomtype adjustment
                    tn = atomtype[n] # do not offset since atomtype is already 1-based
                    sumU[tn+Ne*(k-1)+Ne*K*(m-1)] += c1*c2
                end
            end
        end
    end
end

function threebodydesc!(d3, sumU, pn3, pc3, params)
    nelements = params.counts.nelements
    nabf3 = params.counts.nabf3
    nrbf3 = params.nrbf3
    K3 = params.counts.K3
    Me = nelements*(nelements+1)÷2

    for m in 1:nabf3*nrbf3*Me
        d3[m] = 0.0
    end

    if nelements == 1 #TODO test this 
        for m in 1:nrbf3
            for p in 1:nabf3
                 n1 = pn3[p]+1
                 n2 = pn3[p+1]+1
                 nn = n2 - n1
                for q in 0:nn-1
                    t1 = pc3[n1+q]*sumU[(n1+q) + K3*(m-1)]
                    d3[p + nabf3*(m-1)] += t1*sumU[(n1+q) + K3*(m-1)]
                end
            end
        end
    else
        for m in 1:nrbf3
            for p in 1:nabf3
                n1 = pn3[p]+1
                n2 = pn3[p+1]+1
                nn = n2 - n1
                for q in 0:nn-1
                    k = 0
                    for i1 in 1:nelements
                        t1 = pc3[n1+q]*sumU[i1 + nelements*(n1+q-1) + nelements*K3*(m-1)]
                        for i2 in i1:nelements
                            d3[p + nabf3*(m-1) + nabf3*nrbf3*k] += t1*sumU[i2 + nelements*(n1+q-1) + nelements*K3*(m-1)]
                            k += 1
                        end
                    end
                end
            end 
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

function threebody_forcecoeff!(fb3,cb3,sumU,pn3,pc3,elemindex,params)
    nelements = params.counts.nelements
    nabf3 = params.counts.nabf3
    nrbf3 = params.nrbf3
    K3 = params.counts.K3

    if nelements==1 # TODO Test this case
        for m in 1:nrbf3
            for p in 1:nabf3
                c3 = 2.0 * cb3[p + nabf3*(m-1)]
                n1 = pn3[p]+1
                n2 = pn3[p + 1]+1
                nn = n2 - n1
                idxU = K3*(m-1)
                for q in 0:nn-1
                  k = n1 + q
                  fb3[k + idxU] += c3*pc3[k]*sumU[k + idxU]
                end
            end
        end
    else
        N3 = nabf3 * nrbf3
        for m in 1:nrbf3
            for p in 1:nabf3
                n1 = pn3[p]+1
                n2 = pn3[p + 1]+1
                nn = n2 - n1
                jmp = p + nabf3*(m-1)
                for q in 0:nn-1
                    k = n1 + q  #Combine n1 and q into a single index
                    idxU = nelements*(k-1) + nelements*K3*(m-1)
                    for i1 in 1:nelements
                        tm = pc3[k]*sumU[i1 + idxU]
                        for i2 in i1:nelements
                            em = elemindex[i2+nelements*(i1-1)]
                            t1 = tm*cb3[jmp+N3*em] #em is 0-indexed, jmp has p which is 1-indexed
                            fb3[i2+idxU] += t1
                            fb3[i1+idxU] += pc3[k]*cb3[jmp+N3*em]*sumU[i2+idxU]
                        end
                    end
                end
            end
        end
    end
end

function allbody_forces!(fij,forcecoeff,Ux,Uy,Uz,tj,Nj,params)
    nrbf3 = params.nrbf3
    K3 = params.counts.K3
    nelements = params.counts.nelements
    for j in 1:Nj
        i2 = tj[j] # don't subtract 1, keep 1-indexed
        fx = 0.0
        fy = 0.0
        fz = 0.0
        for m in 1:nrbf3
            for k in 1:K3
                fc = forcecoeff[i2+nelements*(k-1)+nelements*K3*(m-1)]
                idxU = j + Nj*(k-1) + Nj*K3*(m-1)  #Pre-compute the index for abf
                fx += fc * Ux[idxU] 
                fy += fc * Uy[idxU]
                fz += fc * Uz[idxU]
            end
        end
        baseIdx = 3*(j-1)
        fij[baseIdx+1] += fx
        fij[baseIdx+2] += fy
        fij[baseIdx+3] += fz
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
