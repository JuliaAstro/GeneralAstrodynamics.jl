#
# Check Halo.jl
#

using Plots
using Symbolics
using StaticArrays
using LinearAlgebra
using ComponentArrays
using UnitfulAstrodynamics
using DifferentialEquations

wait_for_key(prompt="Press enter to continue...") = (print(stdout, prompt); read(stdin, 1); nothing)

tolerance = reltol = abstol = 1e-14

μ = nondimensionalize(Sun.μ, Earth.μ)
    
r₀ = [0.991975555377273, 0, -0.001871657754011]
v₀ = [0, -0.011750613711503, 0]
τ  = 1.45

Φ  = Matrix{Float64}(undef, 6, 6)

for i ∈ 1:15
    global r₀, v₀, τ

    problem = ODEProblem(
        halo_numerical_tic!,
        ComponentArray(rₛ  = r₀,
                       vₛ  = v₀,
                       Φ₁  = [1.0, 0, 0, 0, 0, 0],
                       Φ₂  = [0, 1.0, 0, 0, 0, 0],
                       Φ₃  = [0, 0, 1.0, 0, 0, 0],
                       Φ₄  = [0, 0, 0, 1.0, 0, 0],
                       Φ₅  = [0, 0, 0, 0, 1.0, 0],
                       Φ₆  = [0, 0, 0, 0, 0, 1.0]),
        (0.0, τ),
        ComponentArray(μ   =  μ)
    )    

    integrator = init(problem, Vern8(); reltol=reltol, abstol=abstol)
    solve!(integrator)

    rₛ = integrator.u.rₛ
    vₛ = integrator.u.vₛ

    Φ = hcat(integrator.u.Φ₁, integrator.u.Φ₂, integrator.u.Φ₃, integrator.u.Φ₄, integrator.u.Φ₅, integrator.u.Φ₆) |> transpose

    MoM = det(Φ)
    evals = eigvals(Φ)

    ∂vₛ = accel(rₛ, vₛ, μ)

    F = @SMatrix [
        Φ[4,3] Φ[4,5] ∂vₛ[1];
        Φ[6,3] Φ[6,5] ∂vₛ[3];
        Φ[2,3] Φ[2,5]  vₛ[2]
    ]

    TERM1 = @SMatrix [r₀[3]; v₀[2]; τ] 
    TERM2 = - inv(F) * @SMatrix [vₛ[1]; vₛ[3]; rₛ[2]] 
    xᵪ = TERM1 + TERM2

    r₀[3] = xᵪ[1]
    v₀[2] = xᵪ[2]
    τ     = xᵪ[3]

    @show F

    @show TERM1
    @show TERM2
    @show xᵪ

    @show Φ

    @show r₀
    @show v₀
    @show τ

    # wait_for_key();
end

