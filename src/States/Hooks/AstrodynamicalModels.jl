#
# Hooks into AstrodynamicalModels.jl
#

"""
$(SIGNATURES)

Returns the `ModelingToolkit.ODESystem` associated 
with `R2BPParameters`, provided by `AstrodynamicalModels`.
"""
model(::R2BPParameters) = AstrodynamicalModels.R2BP

"""
$(SIGNATURES)

Returns the `ModelingToolkit.ODESystem` associated 
with `CR3BPParameters`, provided by `AstrodynamicalModels`.
"""
model(::CR3BPParameters) = AstrodynamicalModels.CR3BP


"""
$(SIGNATURES)

Returns the `DifferentialEquations.ODEFunction` associated 
with `R2BPParameters`, provided by `AstrodynamicalModels`.
"""
vectorfield(::R2BPParameters) = AstrodynamicalModels.R2BPVectorField

"""
$(SIGNATURES)

Returns the `DifferentialEquations.ODEFunction` associated 
with `CR3BPParameters`, provided by `AstrodynamicalModels`.
"""
vectorfield(::CR3BPParameters) = AstrodynamicalModels.CR3BPVectorField