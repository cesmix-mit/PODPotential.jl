using JLD2: load

struct POD_Potential
    basis::PODBasis
    coeff::Vector{Float64}
end

function POD_Potential(basis::PODBasis, coeff_fname::String)
    coeffs = vec(readdlm(coeff_fname, ' '; skipstart=1)) 
    POD_Potential(basis,coeffs)
end

struct LAMMPS_State 
    x::Matrix{Float64}
    type::Vector{Int64}
    map::Vector{Int64}
    inum::Int64
    ilist::Vector{Int64}
    numneigh::Vector{Int64}
    firstneigh::Vector{Vector{Int64}}
    nlocal::Int64
    nghost::Int64
end

function LAMMPS_State(jld_file::String)
    x,type,map,inum,ilist,numneigh,firstneigh,nlocal,nghost = load(jld_file,"x",
                                                              "type",
                                                              "map",
                                                              "inum",
                                                              "ilist",
                                                              "numneigh",
                                                              "firstneigh",
                                                              "nlocal",
                                                              "nghost",)
    return LAMMPS_State(x,type,map,inum,ilist,numneigh,firstneigh, nlocal, nghost)
end


function lammps_compute(state::LAMMPS_State, podpot::POD_Potential)
    
    # state fields are implicitly 0-indexed, I adjust accordingly throughout the code
    x = state.x
    atomtypes = state.type
    inum = state.inum
    ilist = state.ilist
    numneigh = state.numneigh
    firstneigh = state.firstneigh
    elem_map = state.map

    rcut = podpot.basis.params.rcut
    rcutsq = rcut^2

    # pre-compute nijmax    
    nijmax = 0
    for ii in 1:inum
        i = ilist[ii] + 1 # pre-adjust to 1-indexing
        jnum = numneigh[i]

        if nijmax < jnum
            nijmax = jnum
        end
    end

    allocate_workspace!(podpot.basis,nijmax)
    
    fij1 = Vector{Float64}(undef,3*nijmax)
    f    = zeros(Float64,size(x)) 
    
    pot_energy = 0.0 
    forcecoeff = nothing
    for ii in 1:inum
        i = ilist[ii] + 1 # pre-adjust to 1-indexing
        jnum = numneigh[i]

        rij1, ai1, aj1, ti1, tj1, nij = lammpsNeighborList(x,
                                                           firstneigh, 
                                                           atomtypes,
                                                           elem_map,
                                                           numneigh,
                                                           rcutsq,
                                                           i,
                                                           nijmax) # extra argument here
        evdwl = peratomenergyforce2!(fij1,rij1,ti1,tj1,nij,podpot.basis,podpot.coeff)
        pot_energy += evdwl

        tallyforce(f,fij1,ai1,aj1,nij)
    end
    return pot_energy, f
end

function lammpsNeighborList(x, firstneigh, atomtypes,elem_map, numneigh, rcutsq, gi, nijmax)
    # TODO allocation step, need to change to pre-allocate
    ti1 = Vector{Int64}(undef,nijmax)
    tj1 = Vector{Int64}(undef,nijmax)
    ai1 = Vector{Int64}(undef,nijmax)
    aj1 = Vector{Int64}(undef,nijmax)
    rij1 = Vector{Float64}(undef,3*nijmax)
 
    # gi already adjusted to 1-indexing
    nij = 0 # acts as both an index and count, so don't pre-adjust
    itype = elem_map[atomtypes[gi]+1] + 1 # this +1 is in the original code
    ti1[nij+1] = itype; # adjust to 1-indexing
    m = numneigh[gi];
    for l in 1:m
        gj = firstneigh[gi][l] + 1 # pre-adjust to 1-indexing
        delx = x[gj,1] - x[gi,1]
        dely = x[gj,2] - x[gi,2]
        delz = x[gj,3] - x[gi,3]
        rsq = delx*delx + dely*dely + delz*delz
        if rsq < rcutsq && rsq > 1e-20
            rij1[nij*3+1] = delx
            rij1[nij*3+2] = dely
            rij1[nij*3+3] = delz
            ai1[nij+1] = gi 
            aj1[nij+1] = gj 
            ti1[nij+1] = itype
            tj1[nij+1] = elem_map[atomtypes[gj]+1] + 1
            nij += 1
        end
    end 
    return rij1, ai1, aj1, ti1, tj1, nij
end


function tallyforce(force,fij,ai,aj,N)
    for n in 1:N
        im =  ai[n] # these indices already have the 1-index adjustment
        jm =  aj[n] # these indices already have the 1-index adjustment
        nm = 3*(n-1) # adjust for 1-index
        force[im,1] += fij[1+nm]
        force[im,2] += fij[2+nm]
        force[im,3] += fij[3+nm]
        force[jm,1] -= fij[1+nm]
        force[jm,2] -= fij[2+nm]
        force[jm,3] -= fij[3+nm]
    end
end

### Below will be for AtomsCalculator interface...
import AtomsCalculators: potential_energy, 
                        forces, 
                        virial,
                        energy_unit,
                        length_unit



energy_unit(calc::POD_Potential) = u"eV"
length_unit(calc::POD_Potential) = u"Å"

function potential_energy(sys, calc::POD_Potential; 
                          neighbors=nothing,
                          kwargs...) 

end

function forces(sys, calc::POD_Potential;
                neighbors=nothing,
                kwargs...)

end

function virial(sys,calc::POD_Potential;
                neighbors=nothing,
                kwargs...)
    error("virial isn't currently supported for POD_Potential")
end
