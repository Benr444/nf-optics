using Revise
using Images, FileIO, Plots, FFTW
Revise.includet("ROptics.jl")
using Main.ROptics

Plots.default(overwrite_figure = false) # causes display() to make a new window

# the signal begins at x = 0 and ends at x = 10
#signal_length = 10 # meters
function c_signal(x::Real) # x is a length in meters
	k_a = 1 # cycles per meter
	k_b = 2 # cycles per meter
	Φ = 3 * pi / 4
	return sin(2 * pi * k_a * x) + 1//2 * sin(2 * pi * k_b * x + Φ)
end

# signal as a collection of points
Nᵢ = 500 # number of samples for the signal data
kᵢ = 500
xᵢ = 1 / kᵢ
signal_points = [n * xᵢ for n in 0:(Nᵢ - 1)]
signal_data = [c_signal(x) for x in signal_points]

signal_plot = plot(0:(Nᵢ - 1), signal_data, title="signal", label="Nᵢ = $Nᵢ, kᵢ = $kᵢ")
#display(signal_plot)

Nₛ = 8 # number of samples
kₛ = 8 # sampling (spatial) frequency. samples per meter
xₛ = 1 / kₛ # sampling period. meters
sample_points = [n * xₛ for n in 0:(Nₛ - 1)]
sample_data = [c_signal(x) for x in sample_points]

sample_plot = scatter(0:(Nₛ - 1), sample_data, title="sample", label="Nₛ = $Nₛ, kₛ = $kₛ")

# ℱ.signal
ℱ_signal_data = fft(signal_data)
ℱ_signal_plot = plot(0:(Nᵢ - 1), real.(ℱ_signal_data), title="ℱ.signal", label="real")

# ℱ.sample
ℱ_sample_data = fft(sample_data)
ℱ_sample_plot = scatter(
	0:(Nₛ - 1), 
	hcat(
		real.(ℱ_sample_data),
		#imag.(ℱ_sample_data),
	), 
	title="ℱ.sample",
	label=["real" "imag"])

compound_plot = plot(
	signal_plot,
	ℱ_signal_plot,
	sample_plot,
	ℱ_sample_plot,
	#layout = (4, 1),
	framestyle = :origin,
)
display(compound_plot)