using JLD2: load

struct POD_Potential
    basis::PODBasis
    β::Vector{Float64}
end

struct LAMMPS_State 
    x::Matrix{Float64}
    type::Vector{Int64}
    map::Vector{Int64}
    inum::Int64
    ilist::Vector{Int64}
    numneigh::Vector{Int64}
    firstneigh::Vector{Vector{Int64}}
end

function LAMMPS_State(jld_file::String)
   x,type,map,inum,ilist,numneigh,firstneigh = load(jld_file,"x", 
                                                             "type",
                                                             "map",
                                                             "inum",
                                                             "ilist",
                                                             "numneigh",
                                                             "firstneigh")
   return LAMMPS_State(x,type,map,inum,ilist,numneigh,firstneigh)
end


function lammps_compute(state::LAMMPS_State, podpot::POD_Potential)
    
    # everything is implicitly 0-indexed, I adjust accordingly throughout the code
    x = state.x
    atomtypes = state.type
    inum = state.inum
    ilist = state.ilist
    numneigh = state.numneigh
    firstneigh = state.firstneigh
    elem_map = state.map

    rcut = podpot.basis.params.rcut
    rcutsq = rcut^2
    
    nijmax = 0
    for ii in 1:inum
        i = ilist[ii] + 1 # pre-adjust to 1-indexing
        jnum = numneigh[i]

        if nijmax < jnum
            nijmax = jnum
        end

        rij1, ai1, aj1, ti1, tj1 nij, = lammpsNeighborList(x,
                                                           firstneigh, 
                                                           atomtypes,
                                                           elem_map,
                                                           numneigh,
                                                           rcutsq,
                                                           i,
                                                           nijmax) # extra argument here

        evdwl, fij1 = peratomenergyforce2(rij1,ti1,tj1,nij,podpot)

    end
end

function lammpsNeighborList(x, firstneigh, atomtypes,elem_map, numneigh, rcutsq, gi, nijmax)
    # allocation step 
    ti1 = Vector{Int64}(undef,nijmax)
    tj1 = Vector{Int64}(undef,nijmax)
    ai1 = Vector{Int64}(undef,nijmax)
    aj1 = Vector{Int64}(undef,nijmax)
    rij1 = Vector{Float64}(undef,3*nijmax)
 
    # gi already adjusted to 1-indexing
    nij = 0 # acts as both an index and count, so don't pre-adjust
    itype = map[atomtypes[gi]] + 1 # this +1 is in the original code
    ti1[nij+1] = itype; # adjust to 1-indexing
    m = numneigh[gi];
    for l in 1:m
        gj = firstneigh[gi][l] + 1 # pre-adjust to 1-indexing
        delx = x[gj][1] - x[gi][1]
        dely = x[gj][2] - x[gi][2]
        delz = x[gj][3] - x[gi][3]
        rsq = delx*delx + dely*dely + delz*delz
        if rsq < rcutsq && rsq > 1e-20
            rij1[nij*3 +1] = delx
            rij1[nij*3 +2] = dely
            rij1[nij*3 +3] = delz
            ai1[nij+1] = gi 
            aj1[nij+1] = gj 
            ti1[nij+1] = itype
            tj1[nj1+1] = map[atomtypes[gi]] + 1
            nij += 1
        end
    end 
    return rij1, ai1, aj1, ti1, tj1, nij
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
