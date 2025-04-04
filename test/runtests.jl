using PODPotential
using Test
using JLD2 
using BenchmarkTools

@testset "2-body monoclinic HfO2" begin
    pe_ref, f_ref = load("./files/monoclinic_hfo2_2body_output.jld2", "pe", "f")
    
    podbasis = PODBasis([:Hf,:O], 5.5; 
                        rin = 1.0, 
                        besseldegree=4,
                        inversedegree=8,
                        nrbf2=8,
                        nrbf3=0,
                        nrbf4=0,
                        P3=0,
                        P4=0)
    podpot = POD_Potential(podbasis,"./files/simple_2body_HfO2_coefficients.pod")
    lammps_state = LAMMPS_State("./files/monoclinic_hfo2_lammps_state.jld2")
    
    @testset "Correctness" begin 
        pe, f = lammps_compute(lammps_state,podpot)
        
        @test pe ≈ pe_ref
        @test f ≈ f_ref 
    end
end

@testset "3-body monoclinic HfO2" begin
    pe_ref, f_ref = load("./files/monoclinic_hfo2_3body_output.jld2", "pe", "f")
    
    podbasis = PODBasis([:Hf,:O], 5.5; rin = 1.0,
                        besseldegree=4,
                        inversedegree=8,
                        nrbf2=8,
                        nrbf3=6,
                        nrbf4=0,
                        P3=4,
                        P4=0)
    podpot = POD_Potential(podbasis,"./files/simple_3body_HfO2_coefficients.pod")
    lammps_state = LAMMPS_State("./files/monoclinic_hfo2_lammps_state.jld2")
    
    @testset "Correctness" begin 
        pe, f = lammps_compute(lammps_state,podpot)
        
        @test pe ≈ pe_ref
        @test f ≈ f_ref 
    end
end


begin
    println("\n2-body performance:")
    podbasis = PODBasis([:Hf,:O], 5.5; 
                        rin = 1.0, 
                        besseldegree=4,
                        inversedegree=8,
                        nrbf2=8,
                        nrbf3=0,
                        nrbf4=0,
                        P3=0,
                        P4=0)
    podpot = POD_Potential(podbasis,"./files/simple_2body_HfO2_coefficients.pod")
    lammps_state = LAMMPS_State("./files/monoclinic_hfo2_lammps_state.jld2")
    display(@benchmark lammps_compute(lammps_state,podpot))
end

begin
    println("\n3-body performance:")
    podbasis = PODBasis([:Hf,:O], 5.5; rin = 1.0,
                        besseldegree=4,
                        inversedegree=8,
                        nrbf2=8,
                        nrbf3=6,
                        nrbf4=0,
                        P3=4,
                        P4=0)
    podpot = POD_Potential(podbasis,"./files/simple_3body_HfO2_coefficients.pod")
    lammps_state = LAMMPS_State("./files/monoclinic_hfo2_lammps_state.jld2")
    @benchmark lammps_compute(lammps_state,podpot)
end
