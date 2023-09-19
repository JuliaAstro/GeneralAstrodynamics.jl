using OrdinaryDiffEq, AstrodynamicalModels
using SPICE, SPICEBodies, SPICEKernels
using LinearAlgebra
using Plots
using StaticArrays
using Revise, AstrodynamicalSolvers, AstrodynamicalCalculations


plotlyjs()

furnsh(
    de440s(),                   # position and velocity data for major solar system bodies
    latest_leapseconds_lsk(),   # timekeeping, parsing epochs
    gm_de440(),                 # mass parameters for major solar system bodies
    pck00011(),                 # physical properties of major solar system bodies
)

μ = reduced_mass(gm("earth"), gm("moon"))
ic = halo(μ, 1; amplitude=0.01, abstol=1e-14, reltol=1e-14)

u = @MVector [ic.x, 0, ic.z, 0, ic.ẏ, 0]
T = ic.T

halo_ics = let
    problem = ODEProblem(CR3BPFunction(stm=true), MVector{42}(vcat(u, vec(I(6)))), (0, T), (μ,))
    solution = solve(problem, Vern9(), reltol=1e-14, abstol=1e-14, saveat=(T / 100))

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

solution = solve(problem, Vern9(), trajectories=length(manifold_ics), reltol=1e-14, abstol=1e-14)

plot(solution, idxs=(:x, :y, :z), aspect_ratio=1, palette=:rainbow)