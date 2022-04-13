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
#F(x::Matrix{<: Real}) = F(Complex{Float64}.(x))
export F
F⁻¹(x::Matrix{Complex{Float64}}) = ifft(x) # reverse-propagate
#F⁻¹(x::Matrix{<: Real}) = F⁻¹(Complex{Float64}.(x))
export F⁻¹

Plots.png(x::Matrix{Complex{Float64}}, filepath) = png(plot(Gray24.(fit_my_range(real(x)))), filepath)
# ^ method extension to plot
#export png

#set_phase(z, Φ) = complex(abs(z) * cos(Φ), abs(z) * sin(Φ))
set_phase(z, Φ) = abs(z) * exp(1im * Φ)
set_phase(zz::Matrix{Complex{Float64}}, Φs) = broadcast(set_phase, zz, Φs)
export set_phase
#set_modulus(z, l) = complex(abs(l) * cos(angle(z)), abs(l) * sin(angle(z)))
set_modulus(z, r) = abs(r) * exp(1im * angle(z))
set_modulus(zz::Matrix{Complex{Float64}}, rs) = broadcast(set_modulus, zz, rs)
export set_modulus

function gray_modulus(c::AbstractRGBA)
	colors = (red(c), green(c), blue(c))
	(sqrt ∘ sum)((colors.^2))
end
export gray_modulus

# POSSIBLY DEPRECATED?
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

"computes the difference between xx and yy as the sum of squared differences at each value."
function sum_sq_diff(xx, yy)
	sum(broadcast(-, xx, yy).^2)
end
export sum_sq_diff

"rescales every element of xx to fit between 0 and 1"
function cap(xx)
	yy = (xx .- minimum(xx))
	yy ./ maximum(yy)
end
export cap

"computes the mean of xx as the sum of all elements divided by the size of the collection."
function mean(xx)
	sum(xx) / length(xx)
end
export mean

end # MODULE END