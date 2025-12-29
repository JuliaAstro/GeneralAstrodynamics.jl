"""
Tests for orbit types.
"""
module OrbitTests

using AstrodynamicalModels, ModelingToolkit, AstrodynamicalCalculations, Test
using ModelingToolkit: get_p, get_u0

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

    @test system(AstrodynamicalModels.parameters(orbit)) isa ModelingToolkit.System
    @test dynamics(AstrodynamicalModels.parameters(orbit)) isa ModelingToolkit.ODEFunction
    @test system(orbit) isa ModelingToolkit.System
    @test dynamics(orbit) isa ModelingToolkit.ODEFunction
    @test ModelingToolkit.System(orbit) isa ModelingToolkit.System
    @test ModelingToolkit.ODEFunction(orbit) isa ModelingToolkit.ODEFunction
end

end
