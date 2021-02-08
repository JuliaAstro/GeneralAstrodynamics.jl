should_precompile = false

# Don't edit the following! Instead change the script for `snoop_bot`.
ismultios = false
ismultiversion = false
# precompile_enclosure
@static if !should_precompile
    # nothing
elseif !ismultios && !ismultiversion
    @static if (isfile("../deps/SnoopCompile/precompile/precompile_UnitfulAstrodynamics.jl"))
        include("../deps/SnoopCompile/precompile/precompile_UnitfulAstrodynamics.jl")
        _precompile_()
    end
else
    
end # precompile_enclosure
