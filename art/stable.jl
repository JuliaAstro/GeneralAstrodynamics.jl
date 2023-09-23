using GeneralAstrodynamics
using DifferentialEquations, Plots

using Logging
debuglogger = ConsoleLogger(stderr, Logging.Debug)

with_logger(debuglogger) do
    return L1 = halo(EarthMoon, Az = 0.01, L = 1)
end

# L2 = halo(EarthMoon, Az = 0.01, L = 2)

M1 = manifold(
    L1...,
    duration = 1.95 * L1.period,
    eps = -1e-7,
    trajectories = 14,
    direction = Val{:unstable},
)

M2 = manifold(
    L2...,
    duration = 2.1 * L2.period,
    eps = 1e-7,
    trajectories = 14,
    direction = Val{:unstable},
)

artsy = (;
    aspect_ratio = 1.0,
    background = :transparent,
    grid = false,
    axis = nothing,
    xlabel = "",
    ylabel = "",
    zlabel = "",
    border = :none,
    title = "",
    dpi = 600,
)

plot(
    M1;
    vars = :XY,
    artsy...,
    palette = :seaborn_colorblind,
    linestyle = :solid,
    linewidth = 1,
)
plot!(
    M2;
    vars = :XY,
    artsy...,
    palette = :seaborn_colorblind,
    linestyle = :solid,
    linewidth = 1,
)

let xyz = GeneralAstrodynamics.Calculations.secondary_position(massparameter(L1.orbit))
    x, y, z = xyz
    scatter!([x], [y]; markersize = 4, label = :none, color = :black, marker = :x)
end