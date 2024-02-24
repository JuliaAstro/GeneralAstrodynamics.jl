using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ModelingToolkit
using AstrodynamicalModels

x = collect(0.1:0.1:0.6)
p = collect(0.1:0.1:0.7)
t = NaN

model = complete(PlanarEntrySystem())
field = ODEFunction(model)

@show field(x, p, t)

open(joinpath(@__DIR__, "entry.jl"), "w") do file
    func = string(generate_function(model)[1])
    write(file, func)
end

