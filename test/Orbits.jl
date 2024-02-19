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
    p = KeplerianParameters(1e5)

    @test Orbit(u, p) isa Orbit
    @test Orbit(u, p) isa KeplerianOrbit

    u = rand(CartesianState)
    p = R2BParameters(p)

    @test Orbit(u, p) isa Orbit
    @test Orbit(u, p) isa R2BOrbit
end

end