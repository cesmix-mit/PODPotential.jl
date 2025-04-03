using PODPotential
using Test
using JLD2 
using BenchmarkTools

@testset "2-body monoclinic HfO2" begin
    pe_ref, f_ref = load("./files/monoclinic_hfo2_2body_output.jld2", "pe", "f")
    
    podbasis = PODBasis([:Hf,:O], 1.0, 5.5, 4,8,8)
    podpot = POD_Potential(podbasis,"./files/simple_2body_HfO2_coefficients.pod")
    lammps_state = LAMMPS_State("./files/monoclinic_hfo2_lammps_state.jld2")
    
    @testset "Correctness" begin 
        pe, f = lammps_compute(lammps_state,podpot)
        
        @test pe ≈ pe_ref
        @test f ≈ f_ref 
    end

    #@testset "Performance" begin
    #    @benchmark lammps_compute($lammps_state,$podpot)
    #end

end

# Ok, don't know how to get benchmak to show when running test, so doing performance test separately
    podbasis = PODBasis([:Hf,:O], 1.0, 5.5, 4,8,8)
    podpot = POD_Potential(podbasis,"./files/simple_2body_HfO2_coefficients.pod")
    lammps_state = LAMMPS_State("./files/monoclinic_hfo2_lammps_state.jld2")
    @benchmark lammps_compute(lammps_state,podpot)
