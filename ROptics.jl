module ROptics
using Revise
using Images, FileIO, Plots, FFTW
import Base.real

"loads an image file into a matrix of floating-point values."
function make_field(filename)
	image = load(filename)
	return gray_modulus.(image)
end
export make_field

"Fourier-transforms the given matrix." 
F(x::Matrix{Complex{Float64}}) = fft(x) # forward-propagate
#F(x::Matrix{<: Real}) = F(Complex{Float64}.(x))
export F

"Inverse-Fourier-transforms the given matrix." 
F⁻¹(x::Matrix{Complex{Float64}}) = ifft(x) # reverse-propagate
#F⁻¹(x::Matrix{<: Real}) = F⁻¹(Complex{Float64}.(x))
export F⁻¹

"Simplifies the syntax of saving a complex-matrix to file. uses moduli." 
Plots.png(x::Matrix{Complex{Float64}}, filepath) = png(plot(Gray24.(fit_my_range(abs.(x)))), filepath)
# ^ method extension to plot
#export png

"sets the phase of the given complex number z to the given radian-angle Φ."
set_phase(z, Φ) = abs(z) * exp(1im * Φ)
set_phase(zz::Matrix{Complex{Float64}}, Φs) = broadcast(set_phase, zz, Φs)
export set_phase

"sets the modulus of the given complex number z to the given modulus r."
set_modulus(z, r) = abs(r) * exp(1im * angle(z))
set_modulus(zz::Matrix{Complex{Float64}}, rs) = broadcast(set_modulus, zz, rs)
export set_modulus

"converts an image RGB-like value into a float. used in parsing image files."
gray_modulus(c::AbstractRGB) = Float64(Gray(c))
gray_modulus(c::AbstractRGBA) = Float64(Gray(c))
gray_modulus(gc::Gray) = Float64(gc)
export gray_modulus

"deprecated."
function gray_modulus_old(c::AbstractRGBA)
	colors = (red(c), green(c), blue(c))
	(sqrt ∘ sum)((colors.^2))
end
export gray_modulus_old

"#roll the given array so its first position becomes 'centered'."
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
	#Nsq = size(xx)[1] * size(xx)[2]
	sum(abs.(broadcast(-, xx, yy)).^2) / sum(abs.(yy).^2)
end
export f_mse

"rescales every element of xx to fit between 0 and 1"
function cap(xx)
	yy = (xx .- minimum(xx))
	if maximum(yy) == 0.0
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

"selects all elements not in the center of xx. radius is measured in taxicab distance, 
eg radius = 5 will save only (2*5)^2 = 100 elements."
function select_center(xx, radius::Integer)
	yy = copy(xx)
	δx_half = Integer(trunc(size(yy)[1] / 2))
	δy_half = Integer(trunc(size(yy)[1] / 2))
	x_range = (δx_half - radius):(δx_half + radius - 1)
	y_range = (δy_half - radius):(δy_half + radius - 1)
	@show x_range
	@show y_range
	return yy[x_range, y_range]
end
export select_center

"zeros all elements not in the center of xx. radius is measured in taxicab distance, 
eg radius = 5 will save only (2*5)^2 = 100 elements."
function window_center(xx, radius)
	yy = copy(xx)
	δx = size(yy)[1]
	δy = size(yy)[2]
	δx_half = trunc(δx / 2)
	δy_half = trunc(δy / 2)
	for i in 1:δx
		for j in 1:δx
			if abs(i - δx_half) > radius || abs(j - δy_half) > radius
				yy[i, j] = 0
			end
		end
	end
	return yy
end
export select_center

end # MODULE END