#
# Plot recipes for trajectories!
#

@recipe function f(trajectory::Trajectory; vars=:XYZ, timerange=range(first(trajectory.solution.t); stop=last(trajectory.solution.t), length=1000))
    linewidth    --> 1.75
    title        --> "Trajectory" 
    dpi          --> 150
    label        --> :none

    indices = vars isa Tuple ? vars : process_vars(vars)

    labels = Dict(
        0 => "Time" * " ($(timeunit(trajectory)))",
        1 => "X"    * " ($(lengthunit(trajectory)))",
        2 => "Y"    * " ($(lengthunit(trajectory)))",
        3 => "Z"    * " ($(lengthunit(trajectory)))",
        4 => "Ẋ"    * " ($(velocityunit(trajectory)))",
        5 => "Ẏ"    * " ($(velocityunit(trajectory)))",
        6 => "Ż"    * " ($(velocityunit(trajectory)))"
    )

    xguide --> labels[indices[1]] 
    yguide --> labels[indices[2]] 
    if length(indices) == 3
        zguide --> labels[indices[3]]
    end

    getdata(index) = index == 0 ? timerange : map(t -> trajectory(t; idxs=index), timerange)
    
    getdata.((indices...,))
end

@recipe function f(trajectory::AbstractVector{<:Orbit}; vars=:XYZ)
    linewidth    --> 1.75
    title        --> "Trajectory" 
    dpi          --> 150
    label        --> :none

    indices = vars isa Tuple ? vars : process_vars(vars)

    labels = Dict(
        0 => "Time" * " ($(timeunit(first(trajectory))))",
        1 => "X"    * " ($(lengthunit(first(trajectory))))",
        2 => "Y"    * " ($(lengthunit(first(trajectory))))",
        3 => "Z"    * " ($(lengthunit(first(trajectory))))",
        4 => "Ẋ"    * " ($(velocityunit(first(trajectory))))",
        5 => "Ẏ"    * " ($(velocityunit(first(trajectory))))",
        6 => "Ż"    * " ($(velocityunit(first(trajectory))))"
    )

    xguide --> labels[indices[1]] 
    yguide --> labels[indices[2]] 
    if length(indices) == 3
        zguide --> labels[indices[3]]
    end

    getdata(index) = index == 0 ? map(epoch, trajectory) : map(orbit -> state(orbit)[index], trajectory)

    getdata.((indices...,))
end

