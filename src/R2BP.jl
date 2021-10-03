#
# Restricted Two-body Problem models
# 

"""
A `ModelingToolkit.ODESystem` for the Restricted Two-body Problem. 

# Extended Help
The Restricted Two-body Problem is a simplified dynamical model 
describing one small body (spacecraft, etc.) and one celestial 
body. The gravity of the celestial body exhibits a force on the 
small body. This model is commonly used as a simplification to 
descibe our solar systems' planets orbiting our sun, or a 
spacecraft orbiting Earth. 
"""
@memoize function R2BP(; stm=false, structural_simplify=true, name=:R2BP)

    @parameters t μ 
    @variables x(t) y(t) z(t) ẋ(t) ẏ(t) ż(t)
    δ = Differential(t)
    r = @SVector [x,y,z]
    v = @SVector [ẋ,ẏ,ż]

    eqs = vcat(
        δ.(r) .~ v,
        δ.(v) .~ -μ .* (r ./ norm(r)^3)
    )

    if stm 
        @variables Φ[1:6,1:6](t)
        Φ = Symbolics.scalarize(Φ)
        A = Symbolics.jacobian(map(el -> el.rhs, eqs), vcat(r,v))
    
        LHS = map(δ, Φ)
        RHS = map(simplify, A * Φ)

        eqs = vcat(eqs, [LHS[i] ~ RHS[i] for i in 1:length(LHS)])
    end

    if string(name) == "R2BP" && stm 
        modelname = Symbol("R2BPWithSTM")
    else
        modelname = name
    end

    sys = ODESystem(
        eqs, t, stm  ? vcat(r,v,Φ...) : vcat(r,v), [μ]; 
        name = modelname
    )
    return structural_simplify ? ModelingToolkit.structural_simplify(sys) : sys
end

