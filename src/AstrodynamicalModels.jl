"""
Provides astrodynamical models as `ModelingToolkit.ODESystems`. 
Check out the `ModelingToolkit` docs to learn how to use these 
systems for orbit propagation with `DifferentialEquations`, or
see `GeneralAstrodynamics` for some convenient orbit propagation 
wrappers.

# Extended help

## Exports

$(EXPORTS)

## Imports

$(IMPORTS)
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

# Default docstring formatting
using DocStringExtensions
@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(SIGNATURES)

    $(DOCSTRING)
    """

@template (TYPES, CONSTANTS) =
    """
    $(TYPEDEF)

    $(DOCSTRING)
    """

"""
A `ModelingToolkit.ODESystem` for the Restricted Two-body Problem. 


The Restricted Two-body Problem is a simplified dynamical model 
describing one small body (spacecraft, etc.) and one celestial 
body. The gravity of the celestial body exhibits a force on the 
small body. This model is commonly used as a simplification to 
descibe our solar systems' planets orbiting our sun, or a 
spacecraft orbiting Earth. 
"""
const R2BP = let

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

    structural_simplify(R2BP)
end

"""
A `DifferentialEquations`-compatible `ODEFunction` for R2BP dynamics.
Note that this function has several methods, including an in-place 
method! Function signatures follow `ModelingToolkit` and `DifferentialEquations`
conventions. Note that `checkbounds` is __enabled__, so this function 
is not as fast as it could be! Re-call `ODEFunction` with `checkbounds=false`
for an unsafe, but faster, version of this function.
"""
const R2BPVectorField = ODEFunction(R2BP; jac = true, tgrad = false, checkbounds = true, eval_expression=false, eval_module=@__MODULE__)

"""
A `ModelingToolkit.ODESystem` for the Circular Restricted Three-body Problem. 


The Circular Restricted Three-body Problem is a simplified dynamical model 
describing one small body (spacecraft, etc.) and two celestial 
bodies moving in a circle about their common center of mass. 
This may seem like an arbitrary simplification, but this assumption
holds reasonably well for the Earth-Moon, Sun-Earth, and many other 
systems in our solar system.
"""
const CR3BP = let

    @parameters t μ 
    @variables x(t) y(t) z(t) ẋ(t) ẏ(t) ż(t)
    δ = Differential(t)
    r = @SVector [x,y,z]
    v = @SVector [ẋ,ẏ,ż]
              
    eqs = vcat(
        δ.(r) .~ v,
        δ(ẋ)   ~ x + 2ẏ - (μ*(x + μ - 1)*(sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3)) - ((x + μ)*(sqrt(y^2 + z^2 + (x + μ)^2)^-3)*(1 - μ)),
        δ(ẏ)   ~ y - (2ẋ) - (y*(μ*(sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3) + (sqrt(y^2 + z^2 + (x + μ)^2)^-3)*(1 - μ))),
        δ(ż)   ~ z*(-μ*(sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3) - ((sqrt(y^2 + z^2 + (x + μ)^2)^-3)*(1 - μ)))
    )

    @named CR3BP = ODESystem(eqs, t, vcat(r,v), [μ])

    structural_simplify(CR3BP)
end

"""
A `DifferentialEquations`-compatible `ODEFunction` for R2BP dynamics.
Note that this function has several methods, including an in-place 
method! Function signatures follow `ModelingToolkit` and `DifferentialEquations`
conventions. Note that `checkbounds` is __enabled__, so this function 
is not as fast as it could be! Re-call `ODEFunction` with `checkbounds=false`
for an unsafe, but faster, version of this function.
"""
const CR3BPVectorField = ODEFunction(CR3BP; jac = true, tgrad = false, checkbounds = true, eval_expression=false, eval_module=@__MODULE__)

"""
A `ModelingToolkit.ODESystem` for the Circular Restricted Three-body Problem,
with the local linearization included in the state vector and dynamics.
"""
const CR3BPWithSTM = let 

    @parameters t μ 
    @variables x(t) y(t) z(t) ẋ(t) ẏ(t) ż(t) Φ[1:6,1:6](t)
    δ = Differential(t)

    Φ = [
        Φ[i,j] for i in 1:6, j in 1:6 # produces a non-symbolic map
    ]

    r = @SVector [x, y, z]
    v = @SVector [ẋ, ẏ, ż]

    eqs = let

        A = calculate_jacobian(CR3BP)
    
        LHS = map(δ, Φ)
        RHS = map(simplify, A * Φ)

        @assert length(LHS) == length(RHS) == 36 "If this assertion fails, please file an issue at https://github.com/cadojo/AstrodynamicalModels.jl!"

        vcat(
            δ.(r) .~ v,
            δ(ẋ)   ~ x + 2ẏ - (μ*(x + μ - 1)*(sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3)) - ((x + μ)*(sqrt(y^2 + z^2 + (x + μ)^2)^-3)*(1 - μ)),
            δ(ẏ)   ~ y - (2ẋ) - (y*(μ*(sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3) + (sqrt(y^2 + z^2 + (x + μ)^2)^-3)*(1 - μ))),
            δ(ż)   ~ z*(-μ*(sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3) - ((sqrt(y^2 + z^2 + (x + μ)^2)^-3)*(1 - μ))),
            [LHS[i] ~ RHS[i] for i in 1:length(LHS)]
        )

    end

    @named CR3BPWithSTM = ODESystem(eqs, t, vcat(r,v,Φ...), [μ])
    
    structural_simplify(CR3BPWithSTM)
end

"""
A `DifferentialEquations`-compatible `ODEFunction` for R2BP dynamics.
Note that this function has several methods, including an in-place 
method! Function signatures follow `ModelingToolkit` and `DifferentialEquations`
conventions. Note that `checkbounds` is __enabled__, so this function 
is not as fast as it could be! Re-call `ODEFunction` with `checkbounds=false`
for an unsafe, but faster, version of this function.
"""
const CR3BPWithSTMVectorField = ODEFunction(CR3BPWithSTM; jac = false, tgrad = false, checkbounds = true, sparse = false, eval_expression = false, eval_module=@__MODULE__)

end # module
