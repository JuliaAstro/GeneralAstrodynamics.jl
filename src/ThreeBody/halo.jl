#
# Test numerical Halo solver.
# 

using Revise
using UnitfulAstrodynamics, DifferentialEquations, ComponentArrays, LinearAlgebra

μ = nondimensionalize(Sun.μ, Jupiter.μ)
Zₐ = 0.0
ϕ = 0.0u"rad"
L = :L1
hemisphere = :northern
tolerance = 1e-6

r₀, v₀, Τ = halo_analytic(μ; Zₐ=Zₐ, ϕ=ϕ, L=L, hemisphere=hemisphere)
    
μ = 9e-5
tolerance = 1e-6
r₀ = [0.955489,  0.0470234,  0.0]
v₀ = [0.013154, -0.0748174, 0.0]
Τ = 2.9
    
problem = ODEProblem(
    halo_numerical_tic!,
    ComponentArray(rₛ  = [r₀[1], 0, r₀[3]],
                   vₛ  = [0, v₀[2], 0],
                   Φ₁  = copy(I(6)[1,:]),
                   Φ₂  = copy(I(6)[2,:]),
                   Φ₃  = copy(I(6)[3,:]),
                   Φ₄  = copy(I(6)[4,:]),
                   Φ₅  = copy(I(6)[5,:]),
                   Φ₆  = copy(I(6)[6,:])),
    (0.0, Τ),
    ComponentArray(μ   =  μ, 
                   x₁  = -μ, 
                   x₂  =  1 - μ,
                   tol = tolerance)
)

reset = DiscreteCallback(
    (u,t,integrator)->(
        t > 0 && 
        abs(integrator.p.δẋ) ≥ integrator.p.tol &&  
        abs(integrator.p.δż) ≥ integrator.p.tol && 
        u.vₛ[2] == 0.0
    ),
    reset_halo!;
    save_positions=(false, false)
)

sols = solve(problem, callback=reset, reltol=1e-14, abstol=1e-14)