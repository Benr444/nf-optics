module ROptics
using Revise
using Images, FileIO, Plots, FFTW
import Base.real

# ----- OpticField -----
struct OpticParams
	λ::Float64
end
export OpticParams
OpticParams() = OpticParams(1.0)

struct OpticField # the minimum data for an optic-design
	values::Matrix{Complex{Float64}} # 2D Array of Complex numbers
	params::OpticParams
	#scale::Float64 # physical-distance between matrix slots
	#λ::Float64 # physical wavelength
end
export OpticField
OpticField(depths::Matrix{Complex{Float64}}; params = OpticParams()) = OpticField(depths, params)

function OpticField(filename; params = OpticParams())
	image = load(filename)
	gray_array = grayModulus.(image)
	OpticField(gray_array, params)
end

shape(x::OpticField) = shape(x.values)
F(x::OpticField) = OpticField(fft(x.values); x.params) # forward-propagate
export F
F⁻¹(x::OpticField) = OpticField(ifft(x.values), x.params) # reverse-propagate
export F⁻¹

Plots.png(x::OpticField, filepath) = png(plot(Gray24.(fitMyRange(real(x.values)))), filepath)
# ^ method extension to plot
#export png

set_phase(z, Φ) = complex(abs(z) * cos(Φ), abs(z) * sin(Φ))
set_phase(x::OpticField, Φs) = OpticField(broadcast(set_phase, x.values, Φs), x.params)
export set_phase
set_modulus(z, l) = complex(abs(l) * cos(angle(z)), abs(l) * sin(angle(z)))
set_modulus(x::OpticField, Φs) = OpticField(broadcast(set_modulus, x.values, Φs), x.params)
export set_modulus

# ----- END OpticField -----

function grayModulus(c::AbstractRGBA)
	colors = (red(c), green(c), blue(c))
	(sqrt ∘ sum)((colors.^2))
end
export grayModulus

function fitMyRange(xx)
	xx_min = minimum(xx)
	range = maximum(xx) - xx_min
	(xx .- xx_min)./range # -1 -> 1
end
export fitMyRange

end