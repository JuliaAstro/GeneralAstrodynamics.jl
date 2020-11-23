module ThreeBodyUnitTests

using Test
using UnitfulAstrodynamics

@testset "ThreeBody" begin
    
    # Hardcode Gravity parameters for the Sun, 
    # and the Earth-Moon System
    μₛ = 1.32712440018e20u"m^3/s^2"
    μₑ = 4.035032351966808e14u"m^3/s^2"

    # Dimensional initial conditions for spacecraft
    r = [2e9, 7000, 2000]u"km"
    v = [0.001, 0.08, 0.02]u"km/s"
    t = 500u"d"

    # Construct nondimensional state
    sys = ThreeBodySystem(AU, μₛ, μₑ, r, v, t);

    @test true
    
end

