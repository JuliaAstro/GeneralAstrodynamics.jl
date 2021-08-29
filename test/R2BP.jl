"""
Restricted Two-body Model tests.
"""
module R2BPTests 

using LinearAlgebra, Unitful, GeneralAstrodynamics, Test

@testset verbose=false "R2BP Determination" begin
    
    rᵢ = [0.0, 11681.0, 0.0] * u"km"
    vᵢ = [5.134, 4.226, 2.787] * u"km/s"
    corbit = Orbit(CartesianState(rᵢ, vᵢ), Earth)

    @test all(
        keplerian(corbit) .≈ (
            0.723452708202361, 
            24509.265399338536u"km", 
            151.50460766373865u"°", 
            90.0u"°", 
            270.0034742609256u"°", 
            89.99652573907436u"°"
        )
    )
    
    @test CartesianState(cartesian(KeplerianState(keplerian(corbit)...), massparameter(Earth))...) ≈ state(corbit)

    korbit = Orbit(
        KeplerianState(
            0.723452708202361, 
            24509.265399338536u"km", 
            151.50460766373865u"°", 
            90.0u"°", 
            270.0034742609256u"°", 
            89.99652573907436u"°"
        ),
        Earth 
    )

    @test all(
        upreferred.(statevector(KeplerianState(keplerian(cartesian(korbit)..., massparameter(Earth))...))) .≈ 
        upreferred.(statevector(KeplerianState(keplerian(corbit)...)))
    )

end

@testset verbose=false "Kepler's Algorithm" begin

    @testset "Unitful" begin
        
        rᵢ = [0.0, 11681.0, 0.0] * u"km"
        vᵢ = [5.134, 4.226, 2.787] * u"km/s"
        corbit = Orbit(CartesianState(rᵢ, vᵢ), Earth)

        @test all(upreferred.(statevector(state(kepler(corbit)))) .≈ upreferred.(statevector(state(corbit))))

    end

end

@testset verbose=true "Lambert Solvers" begin
    
    @testset "Universal" begin
        
        rᵢ = [0.0, 11681.0, 0.0]u"km"
        vᵢ = [5.134, 4.226, 2.787]u"km/s"
        initial = Orbit(CartesianState(rᵢ, vᵢ), Earth)
    
        Δt = 1000u"s"
        final = kepler(initial, Δt; tol=1e-12)
    
        v₁, v₂ = lambert_universal(position(initial), position(final), massparameter(Earth), Δt; trajectory=:short, tolerance=1e-6, max_iter=1000)

        @test all(v₁ .≈ velocity(initial))
        @test all(v₂ .≈ velocity(final))

    end

    @testset "Oldenhuis" begin
        
        rᵢ = [0.0, 11681.0, 0.0]u"km"
        vᵢ = [5.134, 4.226, 2.787]u"km/s"
        initial = Orbit(CartesianState(rᵢ, vᵢ), Earth)
    
        Δt = 1000u"s"
        final = kepler(initial, Δt; tol=1e-12)
    
        m = 0
        v₁, v₂ = lambert(position(initial), position(final), Δt, m, massparameter(Earth))

        @test_skip all(v₁ .≈ velocity(initial))
        @test_skip all(v₂ .≈ velocity(final))

    end

end


end # module