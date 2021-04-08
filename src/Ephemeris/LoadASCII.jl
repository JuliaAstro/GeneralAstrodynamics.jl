#
# Load ASCII Ephemeris files.
#

"""
NAIF IDs for popular major solar system bodies.
"""
const NAIF_IDS = Dict{String, Int}(
  "solar system barycenter" => 0,
  "sun" => 10,
  "mercury barycenter" => 1,
  "mercury" => 199, 
  "venus barycenter" => 2,
  "venus" => 299,
  "earth-moon barycenter" => 3,
  "moon" => 301,
  "earth" => 399,
  "mars barycenter" => 4,
  "mars" => 499,
  "jupiter barycenter" => 5,
  "jupiter" => 599,
  "saturn barycenter" => 6,
  "saturn" => 699,
  "uranus barycenter" => 7,
  "uranus" => 799,
  "neptune barycenter" => 8,
  "neptune" => 899,
  "pluto barycenter" => 9,
  "pluto" => 999,
)

"""
Loads a comma-delimited ephemeris file and returns a matrix of all values.
"""
function loadascii(filename)
	return hcat(
            map(col->col, 
              filter(col->!all([el=="" for el ∈ col]), 
                collect(eachcol(readdlm(filename, ',')))))...)
end

"""
Given ephemeris file `filename` in CSV format (Julian Day, x, y, z, ẋ, ẏ, ż),
`interpolator` returns a tuple: `tspan`, `rv`.

__Arguments:__
* `data`: the comma-delimited ephemeris data (Cartesian state vectors)
* `type`: the type to convert each element of the data to
* `dateunit`: the units of the first column of the ephemeris data
* `lengthunit`: the units of the 2nd through 4th columns of the ephemeris data
* `velocityunit`: the units of the 5th through 7th columns of the ephemeris data
* `timeunit`: the desired time units of the input to `rv`

__Outputs:__
* `tspan`: a tuple of two Unitful.time quantities, which represent the earliest and 
  latest possible epochs provided by the ephemeris file
* `rv`: a function with one argument (t::Unitful.Time) which returns the body's 
  interpolated Cartesian position and velocity vectors
"""
function interpolator(data; type=Float64, dateunit=u"d", lengthunit=u"km", velocityunit=u"km/s", timeunit=u"s")
  @assert size(data,2) == 7 "Data must be in the following format: [Julian Date, X, Y, Z, Ẋ, Ẏ, Ż]"
	times = convert.(type, ustrip.(timeunit, data[:,1].* dateunit))
	times = range(times[1]; stop=times[end], length=length(times))
	ephem = convert.(type, Matrix(data[:,2:end]))
	
	cart = map(col->CubicSplineInterpolation(times, col), eachcol(ephem))
	
	tspan = (times[1] * timeunit, times[end] * timeunit)
    function rv(t)
        ts = ustrip(timeunit, t)
        
        return (
            SVector{3}(cart[1](ts), cart[2](ts), cart[3](ts)) * lengthunit, 
            SVector{3}(cart[4](ts), cart[5](ts), cart[6](ts)) * velocityunit
        )    
    end
	
    return tspan, rv
end

"""
Wrapper for `interpolator(data; kwargs...)`
"""
interpolator(filename::String; kwargs...) = interpolator(loadascii(filename); kwargs...)