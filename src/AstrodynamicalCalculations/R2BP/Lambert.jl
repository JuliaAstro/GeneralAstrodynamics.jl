#
# Solver for Lambert's problem
#
# References:
# [1] David, A. "Vallado. Fundamentals of Astrodynamics and Applications." (2013).
#

"""
Solves Lambert's problem through the use of univeral variables.
"""
function lambert_universal(r̅₁, r̅₂, μ, Δt; trajectory=:short, tolerance=1e-12, max_iter=25)

    # Specify short way, or long way trajectory
    if trajectory == :short
        tₘ = 1
    elseif trajectory == :long
        tₘ = -1
    else
        throw(ArgumentError("`trajectory` must be set to `:short` or `:long`"))
    end

    r₁ = norm(r̅₁)
    r₂ = norm(r̅₂)

    cosΔν = (r̅₁⋅r̅₂) / (r₁*r₂)
    Δν    = asin(u"rad", tₘ * √(1 - (cosΔν)^2))

    A = upreferred(tₘ * √(r₂*r₁ * (1 + cosΔν)))

    if A ≈ 0
        throw(ErrorException("Can't calculate the orbit."))
    end

    ψₙ = 0.0
    C₂ = 1/2
    C₃ = 1/6

    ψ₊ = 4π^2
    ψ₋ = -4π
    yₙ = r₁ + r₂ + (A * (ψₙ*C₃ - 1) / √(C₂))

    Δtₙ = Δt + 1u"s"
    iter = 0

    while (iter < max_iter) && 
          (abs(Δtₙ - Δt) > (tolerance * oneunit(Δt))) || (A > 0 *oneunit(A) && yₙ < 0 * oneunit(yₙ))

        yₙ = r₁ + r₂ + (A * (ψₙ*C₃ - 1) / √(C₂))
        χₙ = √(yₙ / C₂)
        Δtₙ = (χₙ^3 * C₃ + A*√(yₙ)) / √(μ)

        if Δtₙ < Δt
            ψ₋ = ψₙ
        else
            ψ₊ = ψₙ
        end

        ψₙ = (ψ₊ + ψ₋) / 2
        if ψₙ > tolerance
            C₂ = (1 - cos(√(ψₙ))) /  ψₙ
            C₃ = (√(ψₙ) - sin(√(ψₙ))) / √(ψₙ^3)
        elseif ψₙ < -tolerance
            C₂ = (1 - cosh(√(-ψₙ))) / ψₙ
            C₃ = (sinh(√(-ψₙ)) - √(-ψₙ)) / √((-ψₙ)^3)
        else
            C₂ = 1.0 / 2.0
            C₃ = 1.0 / 6.0
        end

        iter += 1

    end

    f = 1 - yₙ/r₁
    ġ = 1 - yₙ/r₂
    g = A * √(yₙ/μ)

    v̅₁ = upreferred.((r̅₂ .- (f .* r̅₁)) ./ g)
    v̅₂ = upreferred.(((ġ .* r̅₂) .- r̅₁) ./ g)

    return v̅₁, v̅₂
    
end

"""
The following code was converted to Julia, from a [GitHub repository](https://github.com/rodyo/FEX-Lambert)
that hosts a MATLAB implementation. At the time of writing, this respository has a BSD license. I'm providing the copyright notice here, as instructed by the license text.

```
Copyright (c) 2018, Rody Oldenhuis
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of this project.
```
"""
function σmax(y)
    
    an = [
        4.000000000000000e-001,
        2.142857142857143e-001,
        4.629629629629630e-002,
        6.628787878787879e-003,
        7.211538461538461e-004,
        6.365740740740740e-005,
        4.741479925303455e-006,
        3.059406328320802e-007,
        1.742836409255060e-008,
        8.892477331109578e-010,
        4.110111531986532e-011,
        1.736709384841458e-012,
        6.759767240041426e-014,
        2.439123386614026e-015,
        8.203411614538007e-017,
        2.583771576869575e-018,
        7.652331327976716e-020,
        2.138860629743989e-021,
        5.659959451165552e-023,
        1.422104833817366e-024,
        3.401398483272306e-026,
        7.762544304774155e-028,
        1.693916882090479e-029,
        3.541295006766860e-031,
        7.105336187804402e-033
    ]
    
    powers  = y.^(1:25)
    σ       = 4/3 * powers ⋅ an
    ∂σ_∂x   = ((1:25) .* vcat(1, powers[1:24])) ⋅ an
    ∂²σ_∂x² = ((1:25) .* (0:24) .* vcat(1/y, 1, powers[1:23])) ⋅ an
    ∂³σ_∂x³ = ((1:25) .* (0:24) .* (-1:23) .* 
                vcat(1/y/y, 1/y, 1, powers[1:22])) ⋅ an
    
    return σ, ∂σ_∂x, ∂²σ_∂x², ∂³σ_∂x³

end;

"""
The following code was converted to Julia, from a [GitHub repository](https://github.com/rodyo/FEX-Lambert)
that hosts the MATLAB implementation. At the time of writing, the respository has a BSD license. I'm providing the copyright notice here, as instructed by the license text.

```
Copyright (c) 2018, Rody Oldenhuis
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of this project.
```

"""
function LancasterBlanchard(x, q, m)
    
    # Validate input
    if x < -one(x)
        @warn "`x` provided implies negative eccentricity. Correcting..."
        x = abs(x) - 2 * one(x)
    elseif x == -one(x)
        @warn "`x` provided is invalid. Setting `x` to the next float..."
        x = nextfloat(x)
    end
    
    # Parameter E
    E = x*x - one(x)
    
    
    if x == 1 # parabolic input ⟹ analytical solution
        
        T   = 4/3 * (1 - q^3)
        ∂T  = 4/5 * (q^5 - 1)
        ∂²T = ∂T + 120/70 * (1 - q^7)
        ∂³T = 3 * (∂²T - ∂T) + 2400 / 1080 * (q^9 - 1)
        
    elseif abs(x - one(x)) < 1e-2 # almost parabolic ⟹ use series
        
        # Series expansion for T & associated partials
        σ₁, ∂σ_∂x₁, ∂²σ_∂x₂₁, ∂³σ_∂x₃₁ = σmax(-E)
        σ₂, ∂σ_∂x₂, ∂²σ_∂x₂₂, ∂³σ_∂x₃₂ = σmax(-E * q * q)

        T   = σ₁ - q^3 * σ₂
        ∂T  = 2 * x * (q^5 * ∂σ_∂x₂ - ∂σ_∂x₁)
        ∂²T = ∂T/x + 4 * x^2 * (∂²σ_∂x₂₁ - q^7 * ∂²σ_∂x₂₂)
        ∂³T = 3 * (∂²T - ∂T/x) / x + 8 * x * x * (q^9 * ∂³σ_∂x₃₂ - ∂³σ_∂x₃₁)

    else 
        
        y = √(abs(E))
        z = √(1 + q^2 * E)
        f = y * (z - q * x)
        g = x * z - q * E
        
        if E < zero(E)
            δ = atan(f, g) + m*π
        elseif E == zero(E)
            δ = zero(E)
        else
            δ = log(ℯ, max(0, f+g))
        end
        
        T   = 2 * (x - q * z - δ / y) / E

        ∂T  = (4 - 4 * q^3 * x/z - 3 * x * T) / E

        ∂²T = (-4 * q^3/z * (1 - q^2 * x^2 / z^2) - 3*T - 3 * x * ∂T) / E

        ∂³T = (4 * q^3 / z^2 * 
                ((1 - q^2 * x^2 / z^2) + 2 * q^2 * x/z^2 * (z - x)) - 
                    8*∂T - 7*x*∂²T) / E
        
    end
    
    return T, ∂T, ∂²T, ∂³T
    
end

"""
The following code was converted to Julia, from a [GitHub repository](https://github.com/rodyo/FEX-Lambert)
that hosts a MATLAB implementation. At the time of writing, this respository has a BSD license. I'm providing the copyright notice here, as instructed by the license text.

```
Copyright (c) 2018, Rody Oldenhuis
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of this project.
```
"""
function minmax_distances(r̲₁, r̲₂, r₁, r₂, δₜ, a, v̲₁, v̲₂, m, μ)
    
    min_distance = min(r₁, r₂)
    max_distance = max(r₁, r₂)
    
    longway = abs(δₜ) > π # check if longway trajectory was selected
    
    # Find eccentricity vector using triple product identity
    e̲ = ((v̲₁⋅v̲₁) * r̲₁ - (v̲₁ ⋅ r̲₁) * v̲₁) / μ - r̲₁ / r₁
        
    # Scalar eccentricity
    e = norm(e̲)
    
    # Pericenter / apocenter
    pericenter = a * (1 - e)
    apocenter  = Inf # For parabolic / hyperbolic case
    if e < one(e)
        apocenter = a * (1 + e) # elliptical case
    end
    
    # We have the eccentricity vector, so we know exactly
    # where the pericenter is. Given δₜ and that fact 
    # to cross-check if the trajectory goes past it
    
    if m > zero(m) # always elliptical, both apses traversed
        min_distance = pericenter
        max_distance = apocenter
    else # more complicated
        
        pm1 = sign(r₁ * r₁ * (e̲⋅v̲₁) - (r̲₁⋅e̲) * (r̲₁⋅v̲₁))
        pm2 = sign(r₂ * r₂ * (e̲⋅v̲₂) - (r̲₂⋅e̲) * (r̲₂⋅v̲₂))
        
        # Make sure θ₁, θ₂ ∈ (-1, 1)
        θ₁ = pm1 * acos(max(-1, min(1, (r̲₁./r₁) ⋅ (e̲./e))))
        θ₂ = pm2 * acos(max(-1, min(1, (r̲₂./r₂) ⋅ (e̲./e))))
        
        if θ₁+θ₂ < zero(θ₁)
            
            if abs(θ₁) + abs(θ₁) == δₜ   # pericenter was passed
                min_distance = pericenter
            else
                min_distance = apocenter # apocenter was passed
            end
            
        elseif longway
            min_distance = pericenter
            if e < one(e)
                max_distance = apocenter
            end
        end
        
    end
    
    return min_distance, max_distance
    
end

"""
The following code was converted to Julia, from a [GitHub repository](https://github.com/rodyo/FEX-Lambert)
that hosts a MATLAB implementation. At the time of writing, this respository has a BSD license. I'm providing the copyright notice here, as instructed by the license text.

```
Copyright (c) 2018, Rody Oldenhuis
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of this project.
```
"""
function lambert_lancaster_blanchard(
        r̲₁::AbstractVector, r̲₂::AbstractVector, 
        μ::Number, Δt::Number; 
        revolutions = 0,
        branch=:left, 
        trajectory=:short,
        tolerance=1e-12, 
        max_iter=25,
        output_extrema=Val{false})
    
    m = revolutions
    
    if output_extrema == Val{true}
        error_output = [NaN, NaN, NaN] * u"km/s", [NaN, NaN, NaN] * u"km/s", (NaN * u"km", NaN * u"km")
    elseif output_extrema == Val{false}
        error_output = [NaN, NaN, NaN] * u"km/s", [NaN, NaN, NaN] * u"km/s"
    else
        throw(ArgumentError("Keyword argument `output_extrema` must be set to Val{true} or Val{false}."))
    end

    if trajectory ∉ (:short, :long)
        throw(ArgumentError("Must specify :long or :short for trajectory kwarg."))
    end
    
    if Δt < zero(Δt)
        throw(ArgumentError(string(
            "Time of flight must be non-negative!",
            "Use the `trajectory` kwarg to specify",
            "longway or shortway trajectories.")))
    end
    
    if m < zero(m)
        throw(ArgumentError(string(
            "Number of revolutions must be non-negative!",
            "Use the `branch` kwarg to specify",
            "left or righht branch solutions.")))
    end
    
    # Collect unit vectors and magnitudes 
    # of provided position vectors
    r₁ = norm(r̲₁)
    r₂ = norm(r̲₂)
    r̂₁ = normalize(r̲₁)
    r̂₂ = normalize(r̲₂)
    
    # Unit position vector orthogonal to 
    # position vectors
    r̂ₒ = normalize(r̲₁ × r̲₂)

    # Tangential-direction unit vectors
    t̂₁ = r̂ₒ × r̂₁
    t̂₂ = r̂ₒ × r̂₂
    
    # Turn angle 
    δₜ = acos(max(-1, min(1, (r̲₁⋅r̲₂)/r₁/r₂)))
    
    # If longway, account for final velocity
    # direction by flipping δₜ sign
    if trajectory == :long
        δₜ -= 2π
    end
    
    # Constants
    c = √(r₁^2 + r₂^2 - 2*r₁*r₂*cos(δₜ))
    s = (r₁ + r₂ + c) / 2
    T = √(8μ / s^3) * Δt
    q = √(r₁ * r₂) / s * cos(δₜ/2)	
    
    # Initial values
    T₀  = first(LancasterBlanchard(0, q, m))
    δT  = T₀ - T
    phr = mod(2 * atan(1 - q^2, 2 * q), 2π)
    
    # Assume failure
    v₁ = [NaN, NaN, NaN]
    v₂ = v₁
    rₑ = (NaN, NaN)
    
    if m == zero(m) # single revolution
        
        x₀₁ = T₀ * δT / 4/ T
        if δT > zero(Δt)
            x₀ = x₀₁
        else
            x₀₁ = δT / (4 - δT)
            x₀₂ = -√(-δT / (T + T₀/2))
            
            # TODO potential bug, see https://github.com/rodyo/FEX-Lambert/issues/1
            W   = x₀₁ + 1.7 * √(2 - phr/π)  # one of these is right! (1)
            # W   = x₀₁ + 1.7 * √(2 - δₜ/π)     # is it this one?     (2)    
            
            if W ≥ zero(W)
                x₀₃ = x₀₁
            else
                x₀₃ = x₀₁ + ((-W)^(1/16) * (x₀₂ - x₀₁)) # one of these is right! (1)
                # x₀₃ = x₀₁ + (-W^(1/16) * (x₀₂ - x₀₁)) # is it this one?        (2)
            end
            
            λ  = 1 + x₀₃ * (1 + x₀₁)/2 - 0.03 * x₀₃^2 * √(1 + x₀₁)
            x₀ = λ * x₀₃ 
        end
        
        # Check if solution can't be found
        if x₀ < -one(x₀)
            @error "Unable to find solution."
            return error_output
        end
        
    else
        
        # Minumum ∂T(x)
        xMπ = 4 / (3π * (2m + 1))
        if phr < π
            xM₀ = xMπ * (phr / π)^(1/8)
        elseif phr > π
            xM₀ = xMπ * (2 - (2 - phr/π)^(1/8))
        else
            xM₀ = zero(xMπ)
        end
        
        # Halley's method
        xM = xM₀
        ∂T = Inf
        iter = 0
        
        while abs(∂T) > tolerance
            
            iter += 1
            
            _, ∂T, ∂²T, ∂³T = LancasterBlanchard(xM, q, m)
            
            xMp = xM
            xM = xM - 2 * ∂T .* ∂²T ./ (2 * ∂²T.^2 - ∂T .* ∂³T)
            
            if mod(iter, 7) != 0
                xM = (xMp + xM)/2
            end
            
            if iter > max_iter
                @error "Unable to find solution."
                return error_output
            end
            
        end
        
        # xM should be elliptic: xM ∈ (-1, 1)
        if xM < -one(xM) || xM > one(xM)
            @error "Unable to find solution."
            return error_output	
        end
        
        # Corresponding time
        TM = first(LancasterBlanchard(xM, q, m))
        
        # Check that T > minimum
        if TM > T
            @error "Unable to find solution."
            return error_output
        end
        
        # Move onto initial values for second solution!
        
        # Initial values
        TmTM = T - TM
        T0mTM = T₀ - TM
        _, ∂T, ∂²T, __ = LancasterBlanchard(xM, q, m)
        
        # First estimate if m > 0
        if branch == :left
            
            x = √(TmTM / (∂²T/2 + TmTM / (1 - xM)^2))
            W = xM + x
            W = 4 * W / (4 + TmTM) + (1 - W)^2
            x₀  = x * (1 - (1 + m + (δₜ - 1/2)) / 
                    (1 + 0.15m) * x * (W/2 + 0.03x * √W)) + xM
            if x₀ > one(x₀)
                @error "Unable to find solution."
                return error_output
            end
            
        else # Second estimate if m > 0
        
            if δT > zero(δT)
                x₀ = xM - √(TM / (∂²T/2 - TmTM * (∂²T/2/T0mTM - 1/xM^2)))
            else
                x₀₀ = δT / (4 - δT)
                W   = x₀₀ + 1.7 * √Complex(2 * (1 - phr))
                
                if real(W) ≥ zero(real(W))
                    x₀₃ = x₀₀
                else
                    x₀₃ = x₀₀ - √((-W)^(1/8)) * (x₀₀ + √(-δT / (1.5*T₀ - δT)))
                end
                
                W = 4 / (4 - δT)
                λ = (1 + (1 + m + 0.24*(δₜ - 1/2)) / 
                        (1 + 0.15m) * x₀₃ * (W/2 - 0.03x₀₃ * √(W)))
                x₀ = x₀₃ * λ
            end

            if real(x₀) < -one(real(x₀))
                @error "Unable to find solution."
                return error_output
            end
            
        end
        
    end
    
    # Finally, find root of Lancaster & Blancard's function
    
    # Halley's method
    x  = x₀
    Tₓ = Inf
    iter = 0
    
    while abs(Tₓ) > tolerance
        
        iter += 1
        
        Tₓ, ∂T, ∂²T, _ = LancasterBlanchard(x, q, m)
        
        # Find the root of the difference between Tₓ and 
        # required time T
        Tₓ = Tₓ - T
        
        # New value of x
        xₚ = x
        x  = x - 2 * Tₓ * ∂T ./ (2 * ∂T^2 - Tₓ * ∂²T)
        
        if mod(iter, 7) != 0
            x = (xₚ + x) / 2
        end
        
        if iter > max_iter
            @error "Unable to find solution."
            return error_output
        end
        
    end
    
    # Calculate terminal velocities
    
    # Constants
    γ = √(μ * s/2)
    
    if c == zero(c)
        σ = one(x)
        ρ = zero(x)
        z = abs(x)
    else
        σ = 2 * √(r₁*r₂ / c^2) * sin(δₜ/2)
        ρ = (r₁ - r₂) / c
        z = √(1 + q^2 * (x^2 - 1))
    end
    
    # Radial component
    vᵣ₁ =  γ * ((q * z - x) - ρ * (q * z + x)) / r₁
    vᵣ₂ = -γ * ((q * z - x) + ρ * (q * z + x)) / r₂
    v̲ᵣ₁ = vᵣ₁ * r̂₁
    v̲ᵣ₂ = vᵣ₂ * r̂₂
    
    # Tangential component
    vₜ₁ = σ * γ * (z + q*x) / r₁
    vₜ₂ = σ * γ * (z + q*x) / r₂
    v̲ₜ₁  = vₜ₁ * t̂₁
    v̲ₜ₂  = vₜ₂ * t̂₂
    
    # Cartesian velocity
    v̲₁ = v̲ₜ₁ .+ v̲ᵣ₁
    v̲₂ = v̲ₜ₂ .+ v̲ᵣ₂
    
    if output_extrema == Val{true}

        # Find min / max distances
        a = s/2 / (1 - x^2)
        rₘ = minmax_distances(r̲₁, r̲₂, r₁, r₂, δₜ, a, v̲₁, v̲₂, m, μ)
        
        return v̲₁, v̲₂, rₘ

    else

        return v̲₁, v̲₂

    end
    
end


"""
The following code was converted to Julia, from a [GitHub repository](https://github.com/rodyo/FEX-Lambert)
that hosts a MATLAB implementation. At the time of writing, this respository has a BSD license. I'm providing the copyright notice here, as instructed by the license text.

```
Copyright (c) 2018, Rody Oldenhuis
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of this project.
```
"""
function lambert(r1vec, r2vec, tf, m, muC) 

    # original documentation:
    #=
    This routine implements a new algorithm that solves Lambert's problem. The
    algorithm has two major characteristics that makes it favorable to other
    existing ones.

    1) It describes the generic orbit solution of the boundary condition
    problem through the variable X=log(1+cos(alpha/2)). By doing so the
    graph of the time of flight become defined in the entire real axis and
    resembles a straight line. Convergence is granted within few iterations
    for all the possible geometries (except, of course, when the transfer
    angle is zero). When multiple revolutions are considered the variable is
    X=tan(cos(alpha/2)*pi/2).

    2) Once the orbit has been determined in the plane, this routine
    evaluates the velocity vectors at the two points in a way that is not
    singular for the transfer angle approaching to pi (Lagrange coefficient
    based methods are numerically not well suited for this purpose).

    As a result Lambert's problem is solved (with multiple revolutions
    being accounted for) with the same computational effort for all
    possible geometries. The case of near 180 transfers is also solved
    efficiently.

    We note here that even when the transfer angle is exactly equal to pi
    the algorithm does solve the problem in the plane (it finds X), but it
    is not able to evaluate the plane in which the orbit lies. A solution
    to this would be to provide the direction of the plane containing the
    transfer orbit from outside. This has not been implemented in this
    routine since such a direction would depend on which application the
    transfer is going to be used in.

    please report bugs to dario.izzo@esa.int    
    =#

    # adjusted documentation:
    #=
    By default, the short-way solution is computed. The long way solution
    may be requested by giving a negative value to the corresponding
    time-of-flight [tf].

    For problems with |m| > 0, there are generally two solutions. By
    default, the right branch solution will be returned. The left branch
    may be requested by giving a negative value to the corresponding
    number of complete revolutions [m].
    =#

    # Authors
    # .-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.
    # Name       : Dr. Dario Izzo
    # E-mail     : dario.izzo@esa.int
    # Affiliation: ESA / Advanced Concepts Team (ACT)

    # Made more readible and optimized for speed by Rody P.S. Oldenhuis
    # Code available in MGA.M on   http://www.esa.int/gsp/ACT/inf/op/globopt.htm

    # last edited 12/Dec/2009

    # ADJUSTED FOR EML-COMPILATION 24/Dec/2009

    # initial values        
    tol = 1e-14;    bad = false;     days = 86400; 

    # work with non-dimensional units
    r1 = sqrt(dot(r1vec, r1vec));  r1vec = r1vec/r1;
    V = sqrt(muC/r1);          r2vec = r2vec/r1;
    T = r1/V;                  tf    = tf*days/T; # also transform to seconds

    # relevant geometry parameters (non dimensional)
    mr2vec = sqrt(dot(r2vec, r2vec));
    # make 100# sure it's in (-1 <= dth <= +1)
    dth = acos( max(-1, min(1, (dot(r1vec, r2vec))/mr2vec)) );  
        
    # decide whether to use the left or right branch (for multi-revolution
    # problems), and the long- or short way    
    leftbranch = sign(m);   longway = sign(tf);     
    m = abs(m);             tf = abs(tf);
    if (longway < 0)
        dth = 2*pi - dth
    end    

    # derived quantities        
    c      = sqrt(1 + mr2vec^2 - 2*mr2vec*cos(dth)); # non-dimensional chord
    s      = (1 + mr2vec + c)/2;                     # non-dimensional semi-perimeter
    a_min  = s/2;                                    # minimum energy ellipse semi major axis
    Lambda = sqrt(mr2vec)*cos(dth/2)/s;              # lambda parameter (from BATTIN's book)
    crossprd = [r1vec[2]*r2vec[3] - r1vec[3]*r2vec[2],
                r1vec[3]*r2vec[1] - r1vec[1]*r2vec[3], # non-dimensional normal vectors
                r1vec[1]*r2vec[2] - r1vec[2]*r2vec[1]];
    mcr       = sqrt(dot(crossprd, crossprd));           # magnitues thereof
    nrmunit   = crossprd/mcr;                        # unit vector thereof

    # Initial values
    # ---------------------------------------------------------

    # ELMEX requires this variable to be declared OUTSIDE the IF-statement
    logt = log(tf); # avoid re-computing the same value

    # single revolution (1 solution)
    if (m == 0)

        # initial values        
        inn1 = -0.5233;      # first initial guess
        inn2 = +0.5233;      # second initial guess
        x1   = log(1 + inn1);# transformed first initial guess
        x2   = log(1 + inn2);# transformed first second guess

        # multiple revolutions (0, 1 or 2 solutions)
        # the returned soltuion depends on the sign of [m]
    else            
        # select initial values
        if (leftbranch < 0)
            inn1 = -0.5234; # first initial guess, left branch
            inn2 = -0.2234; # second initial guess, left branch
        else
            inn1 = +0.7234; # first initial guess, right branch
            inn2 = +0.5234; # second initial guess, right branch
        end        
        x1 = tan(inn1*pi/2);# transformed first initial guess
        x2 = tan(inn2*pi/2);# transformed first second guess
    end

    # since (inn1, inn2) < 0, initial estimate is always ellipse
    xx   = [inn1, inn2];  aa = a_min ./ (1 .- xx.^2);
    bbeta = longway * 2*asin.(sqrt.((s-c) / 2 ./ aa));
    # make 100.4# sure it's in (-1 <= xx <= +1)
    aalfa = 2*acos.(  max.(-1, min.(1, xx)) );

    # evaluate the time of flight via Lagrange expression
    y12  = @. aa.*sqrt(aa).*((aalfa - sin(aalfa)) - (bbeta-sin(bbeta)) + 2*pi*m);

    # initial estimates for y
    if m == 0
        y1 = log(y12[1]) - logt;
        y2 = log(y12[2]) - logt;
    else
        y1 = y12[1] - tf;
        y2 = y12[2] - tf;
    end

    # Solve for x
    # ---------------------------------------------------------

    # Newton-Raphson iterations
    # NOTE - the number of iterations will go to infinity in case
    # m > 0  and there is no solution. Start the other routine in 
    # that case
    err = Inf;  iterations = 0; xnew = 0;    
    while (err > tol)
        # increment number of iterations
        iterations = iterations + 1;
        # new x
        xnew = (x1*y2 - y1*x2) / (y2-y1);
        # copy-pasted code (for performance)
        if m == 0
            x = exp(xnew) - 1; 
        else 
            x = atan(xnew)*2/pi; 
        end

        a = a_min/(1 - x^2);
        if (x < 1) # ellipse
            beta = longway * 2*asin(sqrt((s-c)/2/a));
            # make 100.4# sure it's in (-1 <= xx <= +1)
            alfa = 2*acos( max(-1, min(1, x)) );
        else # hyperbola
            alfa = 2*acosh(x);
            beta = longway * 2*asinh(sqrt((s-c)/(-2*a)));
        end
        # evaluate the time of flight via Lagrange expression
        if (a > 0)
            tof = a*sqrt(a)*((alfa - sin(alfa)) - (beta-sin(beta)) + 2*pi*m);
        else
            tof = -a*sqrt(-a)*((sinh(alfa) - alfa) - (sinh(beta) - beta));
        end
        # new value of y
        if m ==0
            ynew = log(tof) - logt; 
        else 
            ynew = tof - tf; 
        end

        # save previous and current values for the next iterarion
        # (prevents getting stuck between two values)
        x1 = x2;  x2 = xnew;
        y1 = y2;  y2 = ynew;
        # update error
        err = abs(x1 - xnew);
        # escape clause
        if (iterations > 15)
            bad = true; break; 
        end
    end

    # If the Newton-Raphson scheme failed, try to solve the problem
    # with the other Lambert targeter. 
    if bad
        # NOTE: use the original, UN-normalized quantities
        _branch = leftbranch > 0 ? :left : :right
        _traj   = longway > 0 ? :long : :short
        _m = m

        V1, V2, extremal_distances = lambert_lancaster_blanchard(
            r1vec*r1, r2vec*r1, tf*T, muC;
            revolutions = _m,
            branch = _branch,
            trajectory = _traj,
            output_extrema = Val{true});
        return V1, V2, extremal_distances, 1
    end

    # convert converged value of x
    if m==0
        x = exp(xnew) - 1; 
    else 
        x = atan(xnew)*2/pi; 
    end

    #=
    The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
    now need the conic. As for transfer angles near to pi the Lagrange-
    coefficients technique goes singular (dg approaches a zero/zero that is
    numerically bad) we here use a different technique for those cases. When
    the transfer angle is exactly equal to pi, then the ih unit vector is not
    determined. The remaining equations, though, are still valid.
    =#

    # Solution for the semi-major axis
    a = a_min/(1-x^2);

    # Calculate psi
    if (x < 1) # ellipse
        beta = longway * 2*asin(sqrt((s-c)/2/a));
        # make 100.4# sure it's in (-1 <= xx <= +1)
        alfa = 2*acos( max(-1, min(1, x)) );
        psi  = (alfa-beta)/2;
        eta2 = 2*a*sin(psi)^2/s;
        eta  = sqrt(eta2);
    else       # hyperbola
        beta = longway * 2*asinh(sqrt((c-s)/2/a));
        alfa = 2*acosh(x);
        psi  = (alfa-beta)/2;
        eta2 = -2*a*sinh(psi)^2/s;
        eta  = sqrt(eta2);
    end

    # unit of the normalized normal vector
    ih = longway * nrmunit;

    # unit vector for normalized [r2vec]
    r2n = r2vec/mr2vec;

    # cross-products
    # don't use cross() (emlmex() would try to compile it, and this way it
    # also does not create any additional overhead)
    crsprd1 = [ih[2]*r1vec[3]-ih[3]*r1vec[2],
            ih[3]*r1vec[1]-ih[1]*r1vec[3],
            ih[1]*r1vec[2]-ih[2]*r1vec[1]];    
    crsprd2 = [ih[2]*r2n[3]-ih[3]*r2n[2],
            ih[3]*r2n[1]-ih[1]*r2n[3],
            ih[1]*r2n[2]-ih[2]*r2n[1]];

    # radial and tangential directions for departure velocity
    Vr1 = 1/eta/sqrt(a_min) * (2*Lambda*a_min - Lambda - x*eta);
    Vt1 = sqrt(mr2vec/a_min/eta2 * sin(dth/2)^2);

    # radial and tangential directions for arrival velocity
    Vt2 = Vt1/mr2vec;
    Vr2 = (Vt1 - Vt2)/tan(dth/2) - Vr1;

    # terminal velocities
    V1 = (Vr1*r1vec + Vt1*crsprd1)*V;
    V2 = (Vr2*r2n + Vt2*crsprd2)*V;

    # exitflag
    exitflag = 1; # (success)

    # also compute minimum distance to central body
    # NOTE: use un-transformed vectors again!
    extremal_distances = minmax_distances(r1vec*r1, r2vec*r1, r1, mr2vec*r1, dth, a*r1, V1, V2, m, muC);

    return V1, V2, extremal_distances, exitflag

end

"""
Wrapper for Unitful inputs.
"""
function lambert(r1::AbstractVector{<:Unitful.Length}, r2::AbstractVector{<:Unitful.Length}, tf::Unitful.Time, m::Integer, mu::MassParameter)
    v1, v2, ex, fl = lambert(
        ustrip.(u"km", r1),
        ustrip.(u"km", r2),
        ustrip.(u"d", tf),
        m,
        ustrip.(u"km^3/s^2", mu)
    )

    return v1 .* u"km/s", v2 .* u"km/s", ex .* u"km", fl
end