"""
Tests for CR3BP dynamics.
"""
module CR3BPTests

    using AstrodynamicalModels, ModelingToolkit, Test

    @testset "CR3BP Model Constructors" begin
        model = CR3BP(; stm=false, structural_simplify=false) 
        model = CR3BP(; stm=true, structural_simplify=true) 
        @test true
    end

    @testset "CR3BP Model Calculations" begin
        vectorfield = CR3BPFunction()
        
        r = [0.9253021269565836, 0, 0]
        v = [0.0, 0.05852663414965813, 0.0]
        μ = 0.0009536838895767625
        
        @test isapprox(
                vectorfield(vcat(r,v), [μ], NaN),
                [0.0, 0.05852663414965813,  0.0, 0.053265045303684255, 0.0, -0.0];
                atol=1e-8
        )
    end

end