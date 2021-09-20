"""
Orbit visualization tests. Do plots run without error?
"""
module VisualizationTests

using LinearAlgebra, Unitful, UnitfulAstro, GeneralAstrodynamics, Test
using DifferentialEquations
using Plots

@testset verbose=false "Trajectory Plots" begin
   
    # CR3BP trajectory
    L2Halo = halo(SunEarth; Az = 30_000u"km", L = 2)
    traj   = propagate(L2Halo...)

    try 
        plot(traj)
        @test true
    catch e
        @test false
    end

    # Discrete trajectory
    traj = map(t -> Orbit(traj, t), solution(traj).t)
    try 
        plot(traj)
        @test true
    catch e
        @test false
    end

    # Close all plots
    Plots.closeall()

end

@testset verbose=false "Manifold Plots" begin
   
    # CR3BP manifold
    L2Halo = halo(SunEarth; Az = 30_000u"km", L = 2)

    try 
        plot(manifold(L2Halo...))
        @test true
    catch e
        @test false
    end

    # Close all plots
    Plots.closeall()

end

@testset verbose=false "Energy Plots" begin
   
    # CR3BP trajectory
    L2Halo = halo(EarthMoon; Az = 15_000u"km", L = 2)

    # CR3BP zero velocity curves
    try 
        zerovelocityplot(L2Halo.orbit)
        @test true
    catch e
        @test false
    end

    # Close all plots
    Plots.closeall()

end

end # module