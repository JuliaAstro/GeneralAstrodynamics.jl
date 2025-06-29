module SPICEBodiesExt

using AstrodynamicalModels, SPICEBodies

"""
Return the R2BP parameters associated with the provided body (KM^3/S^2).
"""
AstrodynamicalModels.R2BParameters(body::SPICEBodies.BodyLike) = AstrodynamicalModels.R2BParameters(gm(body))

"""
Return a R2BP orbit with parameters associated with the provided body (KM^3/S^2).
"""
AstrodynamicalModels.R2BOrbit(state, body::SPICEBodies.BodyLike) = AstrodynamicalModels.Orbit(state, R2BParameters(body))

"""
Return the Keplerian parameters associated with the provided body (KM^3/S^2).
"""
AstrodynamicalModels.KeplerianParameters(body::SPICEBodies.BodyLike) = AstrodynamicalModels.KeplerianParameters(gm(body))

"""
Return a Keplerian orbit with parameters associated with the provided body (KM^3/S^2).
"""
AstrodynamicalModels.KeplerianOrbit(state, body::SPICEBodies.BodyLike) = AstrodynamicalModels.Orbit(state, AstrodynamicalModels.KeplerianParameters(body))

"""
Return an orbit with parameters associated with the provided body (KM^3/S^2).
"""
AstrodynamicalModels.Orbit(state, body::SPICEBodies.BodyLike) = state isa OrbitalElements ? AstrodynamicalModels.KeplerianOrbit(state, AstrodynamicalModels.KeplerianParameters(body)) : AstrodynamicalModels.R2BOrbit(state, AstrodynamicalModels.KeplerianParameters(body))


"""
Return a R2BP `ODEProblem` with parameters associated with the provided body (KM^3/S^2).
"""
AstrodynamicalModels.R2BProblem(state, tspan, body::SPICEBodies.BodyLike) = AstrodynamicalModels.R2BProblem(state, tspan, AstrodynamicalModels.R2BParameters(body))

"""
Return the CR3BP parameters associated with the provided body in nondimensional units.
"""
AstrodynamicalModels.CR3BParameters(primary::SPICEBodies.BodyLike, secondary::SPICEBodies.BodyLike) =
    let
        μ₁ = gm(primary)
        μ₂ = gm(secondary)

        μ = min(μ₁, μ₂) / (μ₁ + μ₂)
        AstrodynamicalModels.CR3BParameters(μ)
    end

"""
Return a CR3BP orbit with parameters associated with the provided body in nondimensional units.
"""
AstrodynamicalModels.CR3BOrbit(state, primary::SPICEBodies.BodyLike, secondary::SPICEBodies.BodyLike) = AstrodynamicalModels.Orbit(state, AstrodynamicalModels.CR3BParameters(primary, secondary))

"""
Return a CR3BP `ODEProblem` with parameters associated with the provided bodies in nondimensional units.
"""
AstrodynamicalModels.R2BProblem(state, tspan, primary::SPICEBodies.BodyLike, secondary::SPICEBodies.BodyLike) = AstrodynamicalModels.CR3BProblem(state, tspan, AstrodynamicalModels.CR3BParameters(primary, secondary))

end
