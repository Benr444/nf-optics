module ROptics
using Revise
using Images, FileIO, Plots, FFTW
import Base.real

function make_field(filename)
	image = load(filename)
	return gray_modulus.(image)
end
export make_field

F(x::Matrix{Complex{Float64}}) = fft(x) # forward-propagate
F(x::Matrix{<: Real}) = F(Complex{Float64}.(x))
export F
F⁻¹(x::Matrix{Complex{Float64}}) = ifft(x) # reverse-propagate
F⁻¹(x::Matrix{<: Real}) = F⁻¹(Complex{Float64}.(x))
export F⁻¹

Plots.png(x::Matrix{Complex{Float64}}, filepath) = png(plot(Gray24.(fit_my_range(real(x)))), filepath)
# ^ method extension to plot
#export png

set_phase(z, Φ) = complex(abs(z) * cos(Φ), abs(z) * sin(Φ))
set_phase(x::Matrix{Complex{Float64}}, Φs) = broadcast(set_phase, x, Φs)
export set_phase
set_modulus(z, l) = complex(abs(l) * cos(angle(z)), abs(l) * sin(angle(z)))
set_modulus(x::Matrix{Complex{Float64}}, Φs) = broadcast(set_modulus, x, Φs)
export set_modulus

function gray_modulus(c::AbstractRGBA)
	colors = (red(c), green(c), blue(c))
	(sqrt ∘ sum)((colors.^2))
end
export gray_modulus

function fit_my_range(xx)
	xx_min = minimum(xx)
	range = maximum(xx) - xx_min
	(xx .- xx_min)./range # -1 -> 1
end
export fit_my_range

# roll any array so its first position becomes 'centered'
function half_roll(xx::AbstractArray)
	circshift(xx, size(xx).÷2)
end
export half_roll

end