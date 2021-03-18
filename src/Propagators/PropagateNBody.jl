#
# Propagator for the NBody problem.
#

"""
    nbody_tic

Currently not exported. Used for n-body numerical integration.
"""
function NBodyTic!(∂u, u, p, t=0)
    for i = 1:length(u.body)
        ∂u.body[i].r = u.body[i].v
        ∂u.body[i].v = zero.(∂u.body[i].v)
        for j = 1:length(u.body)
            if i ≠ j
                ∂u.body[i].v += ((p.m[i] * p.m[j]) / norm(u.body[j].r .- u.body[i].r)^3 * (u.body[j].r .- u.body[i].r))
            end    
        end
        ∂u.body[i].v *= (p.G / p.m[i])
    end
    return nothing
end

"""
Uses OrdinaryDiffEq solvers to propagate `sys` Δt into the future.
All keyword arguments are passed directly to OrdinaryDiffEq solvers.

References:
* [1] https://diffeq.sciml.ai/v4.0/tutorials/ode_example.html
* [2] https://github.com/SciML/DifferentialEquations.jl/issues/393#issuecomment-658210231
* [3] https://discourse.julialang.org/t/smart-kwargs-dispatch/14571/15
"""
function propagate(sys::NBodySystem{N,T}, Δt::Unitful.Quantity; kwargs...) where N where T

    # Integration options
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    # Initial conditions
    u₀ = ComponentArray(body=(map(b -> ComponentArray((r=ustrip.(u"m", b.r), v=ustrip.(u"m/s", b.v))), sys.bodies)))
    ts = T.(ustrip.(u"s", (zero(Δt), Δt)))
    p  = ComponentArray((G=6.6743e-11, m=map(b->ustrip(u"kg",b.m), sys.bodies)))

    # Integrate!
    sols = solve(ODEProblem(NBodyTic!, u₀, ts, p); options...)

    # Unpack and return
    bodies = map(x->NBodySystem(
                        map(i->Body(u"m" * x.body[i].r, u"m/s" * x.body[i].v, sys.body[i].m), 
                        1:length(sys.body))), 
                sols.u)

    return Trajectory(
        bodies,
        u"s" * sols.t,
        sols.retcode
    )

end
    
    



    
