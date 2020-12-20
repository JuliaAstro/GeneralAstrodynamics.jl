function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(convert),Type{Float64},CelestialBody{Float64}})
    Base.precompile(Tuple{typeof(convert),Type{Float64},Orbit{Float64}})
    Base.precompile(Tuple{typeof(isapprox),Orbit{Float64},Orbit{Float64}})
    Base.precompile(Tuple{typeof(map),Function,Array{Orbit{Float64},1}})
end
