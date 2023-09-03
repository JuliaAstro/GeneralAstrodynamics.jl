using GeneralAstrodynamics
using DifferentialEquations, Plots

amplitudes = 0.0: 1e-3 : 0.1
halos = map(
    A -> halo(EarthMoon; Az = A, L = 2),
    amplitudes,
)

artsy = (; 
    background = :transparent,
    grid = false,
    axis = nothing,
    xlabel = "", 
    ylabel = "",
    zlabel = "",
    border = :none,
    title = "",
    dpi = 400,
    size = (600, 600),
)

figure = let
    plot(; artsy...)
    for halo in halos
        plot!(
            figure,
            propagate(halo.orbit, 1.1 * halo.period, reltol=1e-14, abstol=1e-14);
            vars = :XYZ,
            artsy...,
            label = :none,
        )
    end
end