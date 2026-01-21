"""
Tests for orbit types.
"""
module OrbitTests

using Test
using AstrodynamicalModels:
    AstrodynamicalModels,
    CartesianState,
    CR3BParameters,
    KeplerianOrbit,
    KeplerianParameters,
    KeplerianState,
    Orbit,
    R2BOrbit,
    R2BParameters,
    dynamics,
    system
using ModelingToolkit: System, ODEFunction, get_p, get_u0

@testset "CartesianState Constructors" begin
    @test rand(CartesianState) isa CartesianState
    @test CartesianState(undef) isa CartesianState
    @test CartesianState() isa CartesianState
    @test CartesianState(randn(3), randn(3)) isa CartesianState
end


@testset "KeplerianState Constructors" begin
    @test rand(KeplerianState) isa KeplerianState
    @test KeplerianState(undef) isa KeplerianState
    @test KeplerianState() isa KeplerianState
end

@testset "Orbit Constructors" begin

    u = rand(KeplerianState)
    p = KeplerianParameters(1e5)

    @test Orbit(u, p) isa Orbit
    @test Orbit(u, p) isa KeplerianOrbit

    u = rand(CartesianState)
    p = R2BParameters(p)

    @test Orbit(u, p) isa Orbit
    @test Orbit(u, p) isa R2BOrbit

end

@testset "ModelingToolkit Integration" begin

    u = rand(CartesianState)
    p = rand(CR3BParameters)
    orbit = Orbit(u, p)

    @test system(AstrodynamicalModels.parameters(orbit)) isa System
    @test dynamics(AstrodynamicalModels.parameters(orbit)) isa ODEFunction
    @test system(orbit) isa System
    @test dynamics(orbit) isa ODEFunction
    @test System(orbit) isa System
    @test ODEFunction(orbit) isa ODEFunction
end

end
