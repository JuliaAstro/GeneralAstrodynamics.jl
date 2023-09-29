module SPICEBodiesExt

using AstrodynamicalModels, SPICEBodies

"""
Return the R2BP parameters associated with the provided body (KM^3/S^2).
"""
AstrodynamicalModels.R2BParameters(body::SPICEBodies.BodyLike) = R2BParameters(gm(body))

"""
Return a R2BP orbit with parameters associated with the provided body (KM^3/S^2).
"""
AstrodynamicalModels.R2BOrbit(state, body::SPICEBodies.BodyLike) = Orbit(state, R2BParameters(body))

"""
Return a R2BP `ODEProblem` with parameters associated with the provided body (KM^3/S^2).
"""
AstrodynamicalModels.R2BProblem(state, tspan, body::SPICEBodies.BodyLike) = R2BProblem(state, tspan, R2BParameters(body))

"""
Return the CR3BP parameters associated with the provided body in nondimensional units.
"""
AstrodynamicalModels.CR3BParameters(primary::SPICEBodies.BodyLike, secondary::SPICEBodies.BodyLike) =
    let
        μ₁ = gm(primary)
        μ₂ = gm(secondary)

        μ = min(μ₁, μ₂) / (μ₁ + μ₂)
        CR3BParameters(μ)
    end

"""
Return a CR3BP orbit with parameters associated with the provided body in nondimensional units.
"""
AstrodynamicalModels.CR3BOrbit(state, primary::SPICEBodies.BodyLike, secondary::SPICEBodies.BodyLike) = Orbit(state, CR3BParameters(primary, secondary))

"""
Return a CR3BP `ODEProblem` with parameters associated with the provided bodies in nondimensional units.
"""
AstrodynamicalModels.R2BProblem(state, tspan, primary::SPICEBodies.BodyLike, secondary::SPICEBodies.BodyLike) = CR3BProblem(state, tspan, CR3BParameters(primary, secondary))

end