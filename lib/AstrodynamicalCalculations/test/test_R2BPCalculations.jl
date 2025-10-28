using AstrodynamicalCalculations:
    AstrodynamicalCalculations,
    argument_of_periapsis,
    cartesian_to_keplerian,
    cartesian_to_perifocal,
    conic,
    inclination,
    kepler,
    keplerian_to_cartesian,
    lambert,
    orbital_period,
    orbital_radius,
    right_ascension_ascending_node,
    specific_energy,
    specific_angular_momentum,
    specific_angular_momentum_vector,
    semimajor_axis

@testset "R2BP Calculations" begin
    @test conic(0) == :Circular
    @test conic(1) == :Parabolic
    @test conic(0.999999) == :Elliptical
    @test conic(1.000001) == :Hyperbolic

    @test collect(cartesian_to_perifocal(
        0.9000000000000004,
        100000.00000000055,
        2.2442063995416452,
        1.1139736859596856,
        5.964076767085505,
        0.786172648028967,
        4.140536323401264,
        -2.892055315559945,
        6.2566039442009,
    )) ≈ collect((
        x = -1.8727989006070282,
        y = 2.645561295501541,
        z = 5.188722983436087,
        ẋ = 7.80864012260457,
        ẏ = -0.8041738032535253,
        ż = 1.7411380870180566,
    ))

    @test semimajor_axis(
        -1.8727989006070282,
        2.645561295501541,
        5.188722983436087,
        7.80864012260457,
        -0.8041738032535248,
        1.7411380870180566,
        398600.4355070226,
    ) == 3.060489599980775

    @test inclination(
        -1.8727989006070282,
        2.645561295501541,
        5.188722983436087,
        7.80864012260457,
        -0.8041738032535248,
        1.7411380870180566,
        398600.4355070226,
    ) == 1.976005846038532

    @test right_ascension_ascending_node(
        -1.8727989006070282,
        2.645561295501541,
        5.188722983436087,
        7.80864012260457,
        -0.8041738032535248,
        1.7411380870180566,
        398600.4355070226,
    ) == 2.943683444332231

    @test argument_of_periapsis(
        -1.8727989006070282,
        2.645561295501541,
        5.188722983436087,
        7.80864012260457,
        -0.8041738032535248,
        1.7411380870180566,
        398600.4355070226,
    ) == 4.316849536012919

    @test specific_angular_momentum_vector(
        -1.8727989006070282,
        2.645561295501541,
        5.188722983436087,
        7.80864012260457,
        -0.8041738032535248,
        1.7411380870180566,
        398600.4355070226,
    ) == specific_angular_momentum_vector(
        -1.8727989006070282,
        2.645561295501541,
        5.188722983436087,
        7.80864012260457,
        -0.8041738032535248,
        1.7411380870180566,
    )

    @test specific_angular_momentum(
        -1.8727989006070282,
        2.645561295501541,
        5.188722983436087,
        7.80864012260457,
        -0.8041738032535248,
        1.7411380870180566,
    ) == norm(
        specific_angular_momentum_vector(
            -1.8727989006070282,
            2.645561295501541,
            5.188722983436087,
            7.80864012260457,
            -0.8041738032535248,
            1.7411380870180566,
        ),
    )

    @test specific_angular_momentum(
        -1.8727989006070282,
        2.645561295501541,
        5.188722983436087,
        7.80864012260457,
        -0.8041738032535248,
        1.7411380870180566,
        398600.4355070226,
    ) == specific_angular_momentum(
        -1.8727989006070282,
        2.645561295501541,
        5.188722983436087,
        7.80864012260457,
        -0.8041738032535248,
        1.7411380870180566,
    )

    @test specific_energy(100000.00000000055, 398600.4355070226) == -1.993002177535102
    @test specific_angular_momentum_vector(
        7501.390893889291,
        7868.393278511291,
        4087.9520744171814,
        4.140536323401264,
        -2.892055315559945,
        6.2566039442009,
        398600.4355070226,
    ) == specific_angular_momentum_vector(
        7501.390893889291,
        7868.393278511291,
        4087.9520744171814,
        4.140536323401264,
        -2.892055315559945,
        6.2566039442009,
    )

    @test specific_angular_momentum(
        7501.390893889291,
        7868.393278511291,
        4087.9520744171814,
        4.140536323401264,
        -2.892055315559945,
        6.2566039442009,
        398600.4355070226,
    ) == norm(
        specific_angular_momentum_vector(
            7501.390893889291,
            7868.393278511291,
            4087.9520744171814,
            4.140536323401264,
            -2.892055315559945,
            6.2566039442009,
        ),
    )

    @test orbital_radius(
        7501.390893889291,
        7868.393278511291,
        4087.9520744171814,
        4.140536323401264,
        -2.892055315559945,
        6.2566039442009,
        398600.4355070226,
    ) == 11279.041186710254

end

@testset "R2BP Determination" begin
    r = [0.0, 11681.0, 0.0]
    v = [5.134, 4.226, 2.787]
    μ = 398600.4354360959

    e, a, i, Ω, ω, ν = cartesian_to_keplerian(r..., v..., μ)
    @test isapprox(
        SVector(e, a, i, Ω, ω, ν),
        SVector(
            0.723452708202361,
            24509.265399338536,
            deg2rad(151.50460766373865),
            π / 2,
            deg2rad(270.0034742609256),
            deg2rad(89.99652573907436),
        ),
        atol = 1e-6,
    )

    x, y, z, ẋ, ẏ, ż = keplerian_to_cartesian(e, a, i, Ω, ω, ν, μ)
    @test isapprox(vcat(r, v), vcat(x, y, z, ẋ, ẏ, ż), atol = 1e-3)
end

@testset "Kepler's Algorithm" begin
    r = [0.0, 11681.0, 0.0]
    v = [5.134, 4.226, 2.787]
    μ = 398600.4354360959
    a = semimajor_axis(r, v, μ)
    T = orbital_period(a, μ)

    (; x, y, z, ẋ, ẏ, ż) = kepler(r..., v..., μ, T)

    @test isapprox(vcat(r, v), vcat(x, y, z, ẋ, ẏ, ż), atol = 1e-3)
end

@testset "Lambert Solvers" begin
    @testset "Universal" begin
        r = [0.0, 11681.0, 0.0]
        v = [5.134, 4.226, 2.787]
        Δt = 1000
        μ = 398600.4354360959

        K = kepler(r..., v..., μ, Δt; atol = 1e-3)

        V1, V2 = lambert(r..., K.x, K.y, K.z, μ, Δt; trajectory = :short, atol = 1e-6)

        @test isapprox(
            vcat(v, K.ẋ, K.ẏ, K.ż),
            vcat(V1.ẋ, V1.ẏ, V1.ż, V2.ẋ, V2.ẏ, V2.ż),
            atol = 1e-3,
        )
    end

    @testset "Lancaster / Blanchard" begin
        r = [0.0, 11681.0, 0.0]
        v = [5.134, 4.226, 2.787]
        Δt = 1000
        μ = 398600.4354360959

        K = kepler(r..., v..., μ, Δt; atol = 1e-12)
        # Unexported
        V1, V2 = AstrodynamicalCalculations.R2BPCalculations.lambert_lancaster_blanchard(
            r,
            SVector(K.x, K.y, K.z),
            Δt,
            μ;
            trajectory = :short,
            atol = 1e-6,
        )

        @test_broken isapprox(
            vcat(v, K.ẋ, K.ẏ, K.ż),
            vcat(V1.ẋ, V1.ẏ, V1.ż, V2.ẋ, V2.ẏ, V2.ż),
            atol = 1e-3,
        )
    end
end
