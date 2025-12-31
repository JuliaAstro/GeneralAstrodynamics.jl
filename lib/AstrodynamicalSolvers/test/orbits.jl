@testset "Orbit Propagation" begin

    orbit = Orbit(rand(CartesianState), rand(R2BParameters))
    @test propagate(orbit, 1e-3) isa ODESolution

    state = copy(orbit.state)
    @test isnothing(propagate!(orbit, 1e-3))
    @test orbit.state != state

end

@testset "Lyapunov Orbit Correction" begin

    μ = 0.012150584395829193
    u = richardson_ic(μ, 1)
    u = AstrodynamicalSolvers.CR3BSolvers.lyapunov(u.x, u.ẏ, μ, u.Δt)


    @test u.x ≈ 0.8222791798525636
    @test u.ẏ ≈ 0.13799313228400178
    @test u.Δt ≈ 2.7536820160579087

end

@testset "Halo Orbit Correction" begin

    μ = 0.012150584395829193
    u = richardson_ic(μ, 2; Z = 0.005)
    u = AstrodynamicalSolvers.CR3BSolvers.halo(u.x, u.z, u.ẏ, μ, u.Δt)

    @test u.x ≈ 1.1202340567932783
    @test u.z ≈ 0.004589679675825104
    @test u.ẏ ≈ 0.17648270824601714
    @test u.Δt ≈ 3.4152029027146815

end

@testset "Semantic Orbit Correction" begin

    μ = 0.012150584395829193
    u = halo(μ, 2)

    @test u.x ≈ 1.124357139749168
    @test u.ẏ ≈ 0.15714566174026515
    @test u.Δt ≈ 3.4068306866985814

    u = halo(μ, 1; amplitude = 0.005)
    @test u.x ≈ 0.823388563881332
    @test u.z ≈ 0.005553604696093592
    @test u.ẏ ≈ 0.12683910108732768
    @test u.Δt ≈ 2.7432058155507653

end

@testset "Monodromy Matrices" begin

    μ = 0.012150584395829193
    ic = halo(μ, 1; amplitude = 0.005)
    u = CartesianState(ic)

    @test monodromy(u, μ, ic.Δt, CR3BFunction(stm = true)) ≈ [
        1317.472125300217 -339.26725145920585 -22.23471304682866 388.2455372345198 126.49047353668445 -3.3779594177227
        -429.7764385389933 111.53280315746214 7.269431052432993 -126.49047353678779 -41.404653982215095 1.0977750840404263
        -11.440384647744368 2.9348175664641527 1.1929568568249131 -3.3779594177162653 -1.0977750840374354 -0.05845992955397675
        3612.12986893901 -929.3271447079647 -60.998106970610365 1064.4911782262984 346.9671305741014 -9.244834479682378
        -1482.5514995781887 382.0700769975138 25.039782109090023 -437.2238230101127 -141.4481439160884 3.8211012689780377
        -75.5369690753493 19.429643984545518 1.2797982252155324 -22.23471304678537 -7.269431052412963 1.1929568568243507
    ]
end
