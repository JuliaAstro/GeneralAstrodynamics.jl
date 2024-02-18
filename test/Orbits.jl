"""
Tests for orbit types.
"""
module OrbitTests

using AstrodynamicalModels, AstrodynamicalCalculations, Test

@testset "CartesianState Constructors" begin
    @test rand(CartesianState) isa CartesianState
    @test CartesianState(undef) isa CartesianState
    @test CartesianState() isa CartesianState
end


@testset "KeplerianState Constructors" begin
    @test rand(KeplerianState) isa KeplerianState
    @test KeplerianState(undef) isa KeplerianState
    @test KeplerianState() isa KeplerianState
end

@testset "Orbit Constructors" begin

    u = rand(KeplerianState)
    p = R2BParameters(1e5)

    @test Orbit(u, p) isa Orbit
    @test Orbit(u, p) isa R2BOrbit
    @test Orbit(u, p) isa KeplerianOrbit

    u = rand(CartesianState)

    @test Orbit(u, p) isa Orbit
    @test Orbit(u, p) isa R2BOrbit
    @test Orbit(u, p) isa CartesianOrbit
end

end