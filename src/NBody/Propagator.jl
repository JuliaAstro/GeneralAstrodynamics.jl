#
# Propagator for the NBody problem.
#

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

function npropagate(sys::System, Δt::Unitful.Quantity, ode_alg::OrdinaryDiffEqAlgorithm = Tsit5(); kwargs...)
   
    # Referencing:
    # [1] https://diffeq.sciml.ai/v4.0/tutorials/ode_example.html
    # [2] https://github.com/SciML/DifferentialEquations.jl/issues/393#issuecomment-658210231
     # [3] https://discourse.julialang.org/t/smart-kwargs-dispatch/14571/15

    # Set default kwargs (modified from [3])
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    u₀ = ComponentArray(body=(map(b -> ComponentArray((r̅=ustrip.(u"m", b.r̅), v̅=ustrip.(u"m/s", b.v̅))), sys.body)))
    tspan = (0.0, ustrip(u"s", Δt))
    params = ComponentArray((G=6.6743e-11, m=map(b->ustrip(u"kg",b.m), sys.body)))

    problem = ODEProblem(nbody_tic, u₀, tspan, params)

    sols = solve(problem, ode_alg; options...)

    body_states = Vector()
    for i = 1:length(u₀.body)
        push!(body_states, PropagationResult(
            u"m" * vcat(map(x->x.body[i].r̅', sols.u)...),
            u"m/s" * vcat(map(x->x.body[i].v̅', sols.u)...)
        ))
    end

    return NBodyPropagationResult(
        u"s" * sols.t,
        body_states,
        sols
    )

end

struct PropagationResult{posType<:Number, velType<:Number}

    r̅::AbstractMatrix{posType}
    v̅::AbstractMatrix{velType}

end

struct NBodyPropagationResult

    t
    body
    ode_solution::ODESolution

end
    
    



    
