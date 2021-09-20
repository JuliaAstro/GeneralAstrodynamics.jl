"""
A module which provides visualizations 
for `Trajectory` and `Manifold` 
instances, as well as other 
common astrodynamical visualizations.

# Extended Help

**Exports**

$(EXPORTS)

**Imports**

$(IMPORTS)
"""
module Visualizations

export zerovelocityplot, zerovelocityplot!

using Plots, RecipesBase
using DifferentialEquations
using ..Calculations
using ..CoordinateFrames
using ..States, ..Propagation
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
Rather than supply indices, users can use `Symbol` 
instances to select indices for plotting. 

# Example

traj = propagate(orbit, period)
plot(traj; vars=:tx)
plot(traj; vars=xyz)
plot(traj; vars=xẋ)
"""
function process_vars(vars::Symbol)
    @assert 2 ≤ length(string(vars)) ≤ 3 "Keyword argument `vars` has an inproper length $(length(vars)). Choose a length in [2, 3]."
    @assert sum(c -> count(c, lowercase(string(vars))), ("t", "x", "y", "z", "ẋ", "ẏ", "ż")) == length(string(vars)) "Invalid character provided in `vars`."

    indices = Dict(
        "t" => 0,
        "x" => 1,
        "y" => 2,
        "z" => 3,
        "ẋ" => 4,
        "ẏ" => 5,
        "ż" => 6
    )

    return [indices[string(char)] for char ∈ lowercase(string(vars))]
end

include(joinpath("Trajectories", "Trajectories.jl"))
include(joinpath("Manifolds", "Manifolds.jl"))
include(joinpath("Energy", "Energy.jl"))

end # module
