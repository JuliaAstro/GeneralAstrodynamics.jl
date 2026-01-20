"""
Tests for R2BP dynamics.
"""
module EntryTests

using Test
using AstrodynamicalModels:
    PlanarEntryFunction,
    PlanarEntryParameters,
    PlanarEntryState,
    PlanarEntrySystem,
    system,
    dynamics
using ModelingToolkit: System, ODEFunction

@testset "Entry Model Constructors" begin
    model = PlanarEntrySystem()
    @test model isa System

    @test rand(PlanarEntryState) isa PlanarEntryState
    @test rand(PlanarEntryParameters) isa PlanarEntryParameters
    @test system(rand(PlanarEntryParameters)) isa System
    @test dynamics(rand(PlanarEntryParameters)) isa ODEFunction

end

@testset "Entry Model Calculations" begin
    vectorfield = PlanarEntryFunction()

    x = [0.1, 0.2, 0.3, 0.4]
    p = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]

    @test isapprox(
        vectorfield(x, p, NaN),
        [1909.7027707373823, -0.7780223808324279, 0.01996668332936563, 0.6633361101853507];
        atol=1e-8
    )
end

end
