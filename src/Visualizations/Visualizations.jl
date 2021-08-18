module Visualizations

using Plots, RecipesBase
using ..Calculations
using ..CoordinateFrames
using ..States, ..Propagation

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

end # module
