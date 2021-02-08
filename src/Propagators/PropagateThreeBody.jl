#
#   PropagateThreeBody.jl
#
#   Includes functions and structures for propagating orbits 
#   within the circular restricted three-body problem.
#

"""
Struct to hold three-body propagation results.
"""
struct ThreeBodyPropagationResult{F<:AbstractFloat} <: PropagationResult 

    t::Vector{F}
    step::Vector{<:NondimensionalThreeBodyState{F}}
    propagation_status::Symbol

    function ThreeBodyPropagationResult(t, step, status)
        T = promote_type(eltype(t), eltype(step[1].rₛ))
        return new{T}(T.(t), [convert(T, step[i]) for i ∈ 1:length(step)], status)
    end

end

"""
Show `ThreebodyPropagationResult` in REPL.
"""
function Base.show(io::IO, result::ThreeBodyPropagationResult)

    println(io, typeof(result), " with ", length(result.t), " timesteps")

end

"""
    threebody_tic

Currently not exported. Used for two-body numerical integration.
"""
function threebody_tic!(∂u, u, p, t)

    ∂u.rₛ =  u.vₛ
    ∂u.vₛ =  accel(u.rₛ, u.vₛ, p.μ)

end

"""
Uses OrdinaryDiffEq solvers to propagate `orbit` Δt into the future.
All keyword arguments are passed directly to OrdinaryDiffEq solvers.
"""
function propagate(sys::NondimensionalThreeBodyState, 
                   Δt, 
                   ode_alg::OrdinaryDiffEqAlgorithm=Tsit5();
                   kwargs...)

    # Referencing:
    # [1] https://diffeq.sciml.ai/v4.0/tutorials/ode_example.html
    # [2] https://github.com/SciML/DifferentialEquations.jl/issues/393#issuecomment-658210231
    # [3] https://discourse.julialang.org/t/smart-kwargs-dispatch/14571/15

    # Set default kwargs (modified from [3])
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    # Define the problem (modified from [2])
    problem = ODEProblem(threebody_tic!, 
                         ComponentArray((rₛ=sys.rₛ, vₛ=sys.vₛ)), 
                         (0.0, Δt), 
                         ComponentArray((μ=sys.μ, x₁=sys.r₁[1], x₂=sys.r₂[1])))

    # Solve the problem! 
    sols = solve(problem, ode_alg; options...)

    # Return PropagationResult structure
    return ThreeBodyPropagationResult(
        sols.t,
        map(step->NondimensionalThreeBodyState(step.rₛ, step.vₛ, sys.μ, sys.DU, sys.DT), sols.u),
        sols.retcode
    )

end

