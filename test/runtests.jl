"""
Restricted Two-body Model tests.
"""
module AstrodynamicalCalculationsTests

using Test, AstrodynamicalCalculations
using StaticArrays

@testset verbose=false "R2BP Determination" begin
    
    r = [0.0, 11681.0, 0.0] 
    v = [5.134, 4.226, 2.787] 
    μ = 398600.4354360959

    e, a, i, Ω, ω, ν = cartesian_to_keplerian(r, v, μ)
    @test isapprox(
        SVector(e, a, i, Ω, ω, ν), 
        SVector(
            0.723452708202361, 
            24509.265399338536, 
            deg2rad(151.50460766373865), 
            π/2, 
            deg2rad(270.0034742609256), 
            deg2rad(89.99652573907436)
        ), 
        atol = 1e-6,
    )
    
    rₙ, vₙ = keplerian_to_cartesian(e, a, i, Ω, ω, ν, μ)
    @test isapprox(vcat(r, v),  vcat(rₙ, vₙ), atol=1e-3)

end

@testset verbose=false "Kepler's Algorithm" begin
        
    r = [0.0, 11681.0, 0.0] 
    v = [5.134, 4.226, 2.787]
    μ = 398600.4354360959
    a = semimajor_axis(r, v, μ)
    T = orbital_period(a, μ)

    rₙ, vₙ = kepler(r, v, μ, T)

    @test isapprox(vcat(r,v), vcat(rₙ, vₙ), atol=1e-3)
end

@testset verbose=true "Lambert Solvers" begin
    
    @testset "Universal" begin
        
        r = [0.0, 11681.0, 0.0]
        v = [5.134, 4.226, 2.787]
        Δt = 1000
        μ = 398600.4354360959

        rₙ, vₙ = kepler(r, v, Δt; atol=1e-3)
    
        v₁, v₂ = lambert(r, rₙ, μ, Δt; trajectory=:short, atol=1e-6)

        @test isapprox(vcat(v, vₙ), vcat(v₁, v₂), atol=1e-3)

    end

    @testset "Oldenhuis" begin
        r = [0.0, 11681.0, 0.0]
        v = [5.134, 4.226, 2.787]
        Δt = 1000
        μ = 398600.4354360959

        rₙ, vₙ = kepler(r, v, Δt; atol=1e-12)
    
        v₁, v₂ = lambert_oldenhuis(r, rₙ, μ, Δt; trajectory=:short, atol=1e-6)

        @test_broken isapprox(vcat(v, vₙ), vcat(v₁, v₂), atol=1e-3)
    end

end

@testset "CR3BP Calculations" begin
    r = [ 
        1.007988
        0.0
        0.001864
    ]

    v = zeros(3,1)

    μ = 3.003480593992993e-6

    @test jacobi_constant(r, v, μ) ≈ 3.000907212196274
end

end # module