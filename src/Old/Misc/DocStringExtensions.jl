# 
# Provides common DocStringExtensions abbreviations
# for UnitfulAstrodynamics.jl sub-modules. 
#

@template (FUNCTIONS, METHODS, MACROS) =
"""
$(SIGNATURES)
$(DOCSTRING)
"""

@template (DEFAULT) = 
"""
$(DOCSTRING)
"""