"""
Tests for R2BP dynamics.
"""
module EntryTests

using AstrodynamicalModels, ModelingToolkit, Test
using ModelingToolkit: get_p, get_u0

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
    sys = vectorfield.sys

    u0 = get_u0(sys, [
        :γ => 0.1
        :v => 0.2
        :r => 0.3
        :θ => 0.4
    ])

    p = get_p(sys, [
        :R => 0.1
        :P => 0.2
        :H => 0.3
        :m => 0.4
        :A => 0.5
        :C => 0.6
        :μ => 0.7
    ])

    @test isapprox(
        vectorfield(u0, p, NaN),
        [1909.7027707373823, -0.7780223808324279, 0.01996668332936563, 0.6633361101853507];
        atol=1e-8
    )
end

end
