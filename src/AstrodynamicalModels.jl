"""
Provides astrodynamical models as `ModelingToolkit.ODESystems`. 
Check out the `ModelingToolkit` docs to learn how to use these 
systems for orbit propagation with `DifferentialEquations`, or
see `GeneralAstrodynamics` for some convenient orbit propagation 
wrappers.
"""
module AstrodynamicalModels

# NOTE right at the top. DO NOT CHANGE THE ORDER
# OF THE VARIABLES IN THE MODELS BELOW. 
# Downstream users do not have any way to 
# specify individual states when constructing 
# an `ODEProblem` from each `ODESystem`.

# Export every model!
export R2BP, CR3BP, CR3BPWithSTM
export R2BPVectorField, CR3BPVectorField, CR3BPWithSTMVectorField

# AstrodynamicalSystems.jl simply defines 
# `*System` variables that represent common 
# astrodynamical models: R2BPP, CR3BP, etc.
using Symbolics, ModelingToolkit

# Provides the `norm` function
using LinearAlgebra

# Provides @SVector
using StaticArrays

# Allows for providing compiled functions to user
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)


"""
A `ModelingToolkit.ODESystem` for the Restricted Two-body Problem. 


The Restricted Two-body Problem is a simplified dynamical model 
describing one small body (spacecraft, etc.) and one celestial 
body. The gravity of the celestial body exhibits a force on the 
small body. This model is commonly used as a simplification to 
descibe our solar systems' planets orbiting our sun, or a 
spacecraft orbiting Earth. 
"""
R2BP = let

    @parameters t μ 
    @variables x(t) y(t) z(t) ẋ(t) ẏ(t) ż(t)
    δ = Differential(t)
    r = @SVector [x,y,z]
    v = @SVector [ẋ,ẏ,ż]

    eqs = vcat(
        δ.(r) .~ v,
        δ.(v) .~ -μ .* (r ./ norm(r)^3)
    )

    @named R2BP = ODESystem(eqs, t, vcat(r,v), [μ])

end

"""
A `DifferentialEquations`-compatible `ODEFunction` for R2BP dynamics.
Note that this function has several methods, including an in-place 
method! Function signatures follow `ModelingToolkit` and `DifferentialEquations`
conventions.
"""
const R2BPVectorField = ODEFunction(R2BP; jac = true, tgrad = true, eval_expression=false, eval_module=@__MODULE__)

"""
A `ModelingToolkit.ODESystem` for the Circular Restricted Three-body Problem. 


The Circular Restricted Three-body Problem is a simplified dynamical model 
describing one small body (spacecraft, etc.) and two celestial 
bodies moving in a circle about their common center of mass. 
This may seem like an arbitrary simplification, but this assumption
holds reasonably well for the Earth-Moon, Sun-Earth, and many other 
systems in our solar system.
"""
CR3BP = let

    @parameters t μ 
    @variables x(t) y(t) z(t) ẋ(t) ẏ(t) ż(t)
    δ = Differential(t)
    r = @SVector [x,y,z]
    v = @SVector [ẋ,ẏ,ż]

    eqs = vcat(
        δ.(r) .~ v,
        δ(ẋ) ~ x + 2ẏ - (μ*(μ + x - 1)*(sqrt((μ + x - 1)^2 + y^2 + z^2)^-3)) - ((μ + x)*(sqrt(y^2 + z^2 + (μ + x)^2)^-3)*(1 - μ)),
        δ(ẏ) ~ y - (ẋ) - (y*(μ*(sqrt((μ + x - 1)^2 + y^2 + z^2)^-3) + (sqrt(y^2 + z^2 + (μ + x)^2)^-3)*(1 - μ))),
        δ(ż) ~ z*(-μ*(sqrt((μ + x - 1)^2 + y^2 + z^2)^-3) - ((sqrt(y^2 + z^2 + (μ + x)^2)^-3)*(1 - μ)))
    )

    @named CR3BP = ODESystem(eqs, t, vcat(r,v), [μ])

end

"""
A `DifferentialEquations`-compatible `ODEFunction` for R2BP dynamics.
Note that this function has several methods, including an in-place 
method! Function signatures follow `ModelingToolkit` and `DifferentialEquations`
conventions.
"""
const CR3BPVectorField = ODEFunction(CR3BP; jac = true, tgrad = true, eval_expression=false, eval_module=@__MODULE__)

"""
A `ModelingToolkit.ODESystem` for the Circular Restricted Three-body Problem,
with the local linearization included in the state vector and dynamics.
"""
CR3BPWithSTM = let 

    @parameters t μ 
    @variables Φ[1:6,1:6](t)
    δ = Differential(t)
    
    x, y, z, ẋ, ẏ, ż = states(CR3BP)

    r = @SVector [x, y, z]
    v = @SVector [ẋ, ẏ, ż]
    
    # Potential energy of spacecraft (this code was symbolically generated with Symbolics.jl, then copied and pasted)
    U = μ*(sqrt((μ + x - 1)^2 + y^2 + z^2)^-1) + (1//2) * (x^2) + (1//2) * (y^2) + (sqrt(y^2 + z^2 + (μ + x)^2)^-1)*(1 - μ)

    # Hessian of potential energy
    H = Symbolics.hessian(U, r)

    eqs = let
        A = vcat(
            Matrix{Float64}(hcat(zeros(3,3), I(3))),
            hcat(H, Float64[0 2 0; -2 0 0; 0 0 0])
        )

        LHS = [
            δ(Φ[i,j]) for j in 1:6 for i in 1:6
        ]

        RHS = [
            A[i,j] * Φ[i,j] for j in 1:6 for i in 1:6
        ]

        @assert length(LHS) == length(RHS) == 36 "If this assertion fails, please file an issue at https://github.com/cadojo/AstrodynamicalModels.jl!"

        [LHS[i] ~ RHS[i] for i in 1:length(LHS)]
    end

    @named CR3BPWithSTM = ODESystem(vcat(equations(CR3BP), eqs), t, vcat(r,v,[col for col in eachcol(Φ)]...), [μ])
end

"""
A `DifferentialEquations`-compatible `ODEFunction` for R2BP dynamics.
Note that this function has several methods, including an in-place 
method! Function signatures follow `ModelingToolkit` and `DifferentialEquations`
conventions.
"""
const CR3BPWithSTMVectorField = ODEFunction(CR3BPWithSTM; jac = true, tgrad = false, sparse = false, eval_expression = false, eval_module=@__MODULE__)

end # module
