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
  df = DataFrame(CSV.File(filename))
  return hcat(map(col->col, filter(col->!all(x->ismissing(x), col), collect(eachcol(df))))...)
end

"""
A structure for storing a `CubicSplineInterpolation`
for an ephemeris file.
"""
struct Interpolator{T1<:Unitful.Time, T2<:Unitful.Time, F<:Function}
  timespan::Tuple{T1,T2}
  state::F
  Interpolator(tspan, func) = new{typeof(tspan[1]), typeof(tspan[2]), typeof(func)}(tspan, func)
end

"""
Prints an `Interpolator` instance to `io`.
"""
Base.show(io::IO, interp::Interpolator) = println(io, "Interpolated `CartesianState` ephemeris data. Valid within the following Julian epochs: [$(string(interp.timespan[1])), $(string(interp.timespan[2]))].")

"""
Prints an `Interpolator` instance to `io`.
"""
Base.show(io::IO, ::MIME"text/plain", interp::Interpolator) = show(io, interp)

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
  function rv(t; frame=Inertial)
      ts = ustrip(timeunit, t)
      
      return CartesianState(
          SVector{3}(cart[1](ts), cart[2](ts), cart[3](ts)) * lengthunit, 
          SVector{3}(cart[4](ts), cart[5](ts), cart[6](ts)) * velocityunit,
          t, Inertial
      )    
  end

  return Interpolator(tspan, (t; frame=Inertial) -> let
        ts = ustrip(timeunit, t)
        CartesianState(
            SVector{3}(cart[1](ts), cart[2](ts), cart[3](ts)) * lengthunit, 
            SVector{3}(cart[4](ts), cart[5](ts), cart[6](ts)) * velocityunit,
            t, Inertial
        )    
    end
  )
end

"""
Wrapper for `interpolator(data; kwargs...)`
"""
interpolator(filename::String; kwargs...) = interpolator(loadascii(filename); kwargs...)