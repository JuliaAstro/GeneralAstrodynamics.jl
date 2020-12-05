#
# Propagator for the NBody problem.
#

"""
Struct to hold n-body propagation results.
"""
struct MultibodyPropagationResult{T,B} <: PropagationResult

    t::T
    step::B
    ode_solution::ODESolution

end

"""
Show `MultibodyPropagationResult` in REPL.
"""
function Base.show(io::IO, result::MultibodyPropagationResult)

    println(io, typeof(result), " with ", length(result.t), " timesteps")
    println(io, "  ", "t::", string(typeof(result.t)))
    println(io, "  ", "step::", string(typeof(result.step)))
    println(io, "  ", "ode_solution::", string(typeof(result.ode_solution)))

end

"""
    nbody_tic

Currently not exported. Used for n-body numerical integration.
"""
function nbody_tic(u, p, t)

    ∂u = ComponentArray(body = map(x->ComponentVector(r̅=SVector(0.0, 0.0, 0.0), v̅=SVector(0.0, 0.0, 0.0)), u.body))
    for i = 1:length(u.body)

        ∂u.body[i].r̅ = u.body[i].v̅
        for j = 1:length(u.body)

            if i ≠ j
                ∂u.body[i].v̅ += ((p.m[i] * p.m[j]) / norm(u.body[j].r̅ .- u.body[i].r̅)^3 * (u.body[j].r̅ .- u.body[i].r̅))
            end    

        end
        ∂u.body[i].v̅ *= (p.G / p.m[i])

    end

    return ∂u

end

"""
Uses OrdinaryDiffEq solvers to propagate `sys` Δt into the future.
All keyword arguments are passed directly to OrdinaryDiffEq solvers.
"""
function propagate(sys::MultibodySystem, Δt::Unitful.Quantity, ode_alg::OrdinaryDiffEqAlgorithm = Tsit5(); kwargs...)
   
    if !all(map(body->body.m ≥ 0u"kg", sys.body))
        @warn "One or more bodies have a mass that is not ≥ 0: this is not supported, and the propagation result will not be correct."
    end
    
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    u₀ = ComponentArray(body=(map(b -> ComponentArray((r̅=ustrip.(u"m", b.r̅), v̅=ustrip.(u"m/s", b.v̅))), sys.body)))
    tspan = (0.0, ustrip(u"s", Δt))
    params = ComponentArray((G=6.6743e-11, m=map(b->ustrip(u"kg",b.m), sys.body)))

    problem = ODEProblem(nbody_tic, u₀, tspan, params)
    sols = solve(problem, ode_alg; options...)

    bodies = map(x->MultibodySystem(
                        map(i->Body(u"m" * x.body[i].r̅, u"m/s" * x.body[i].v̅, sys.body[i].m), 
                        1:length(sys.body))), 
                    sols.u)

    return MultibodyPropagationResult(
        u"s" * sols.t,
        bodies,
        sols
    )

end
    
    



    
