using OrdinaryDiffEq, AstrodynamicalModels
using SPICE, SPICEBodies, SPICEKernels
using AstrodynamicalCalculations
using LinearAlgebra
using Plots
using Revise, AstrodynamicalSolvers

plotlyjs()

furnsh(
    de440s(),                   # position and velocity data for major solar system bodies
    latest_leapseconds_lsk(),   # timekeeping, parsing epochs
    gm_de440(),                 # mass parameters for major solar system bodies
    pck00011(),                 # physical properties of major solar system bodies
)

μ = reduced_mass(gm("earth"), gm("moon"))

using Logging
debuglogger = ConsoleLogger(stderr, Logging.Debug)

with_logger(debuglogger) do
    u, T = halo(μ, 1; amplitude=0.01)
end

halo_ics = let
    problem = ODEProblem(CR3BPFunction(stm=true), vcat(u, vec(I(6))), (0, T), (μ,))
    solution = solve(problem, Vern9(), reltol=1e-12, abstol=1e-12, saveat=T / 20)

    solution.u
end

Φ = monodromy(u, μ, T)

manifold_ics = [
    diverge(u[begin:begin+5], reshape(u[begin+6:end], 6, 6), Φ; eps=-1e-7)
    for u in halo_ics
]

problem = EnsembleProblem(
    ODEProblem(CR3BPFunction(), u, (0.0, 1.95T), (μ,)),
    prob_func=(prob, i, repeat) -> remake(prob; u0=manifold_ics[i]),
)

solution = solve(problem, Vern9(), trajectories=14, reltol=1e-12, abstol=1e-12)

plot(solution, idxs=(:x, :y))