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

set_phase(z, Φ) = abs(z) * exp(1im * Φ)
set_phase(zz::Matrix{Complex{Float64}}, Φs) = broadcast(set_phase, zz, Φs)
export set_phase
set_modulus(z, r) = abs(r) * exp(1im * angle(z))
set_modulus(zz::Matrix{Complex{Float64}}, rs) = broadcast(set_modulus, zz, rs)
export set_modulus

gray_modulus(c::AbstractRGBA) = Float64(Gray(c))
export gray_modulus
function gray_modulus_old(c::AbstractRGBA)
	colors = (red(c), green(c), blue(c))
	(sqrt ∘ sum)((colors.^2))
end
export gray_modulus_old

# roll any array so its first position becomes 'centered'
function half_roll(xx::AbstractArray)
	circshift(xx, size(xx).÷2)
end
export half_roll

"computes the difference between xx and yy as the sum of squared differences at each value."
function sum_sq_diff(xx, yy)
	sum(abs.(broadcast(-, xx, yy)).^2)
end
export sum_sq_diff

"computes the mean squared error between xx and yy as the average sum of squared differences at each value."
function mean_se(xx, yy)
	mean(abs.(broadcast(-, xx, yy)).^2)
end
export mean_se

"computes the 2-norm (euclidean) between xx and yy as the sqrt of sum of squared differences at each value."
function norm(xx, yy)
	sqrt(sum_sq_diff(xx, yy))
end
norm(xx) = norm(xx, 0) # non-difference norm.
export norm

"fienup-defined mean-squared-error of xx and yy."
function f_mse(xx, yy)
	sum(abs.(broadcast(-, xx, yy)).^2) / sum(abs.(yy).^2)
end
export f_mse

"rescales every element of xx to fit between 0 and 1"
function cap(xx)
	yy = (xx .- minimum(xx))
	if maximum(yy) == 0.0
		@show maximum(yy)
		return yy
	else
		return yy ./ maximum(yy)
	end
end
export cap

"computes the mean of xx as the sum of all elements divided by the size of the collection."
function mean(xx)
	sum(xx) / length(xx)
end
export mean

"square: rounds x to 3 sigfigs. useful for pipe notation 3.999... |> □."
□(x) = round(x, sigdigits = 3)
export □

"blacksquare: returns a rounding function to n significant figures. default is 3."
■(n = 3) = (x -> round(x, sigdigits = n))
export ■

end # MODULE END