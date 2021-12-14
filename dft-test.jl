using Revise
using Images, FileIO, Plots, FFTW
Revise.includet("ROptics.jl")
using Main.ROptics

# see page 54 of lyon's "understanding digital signal processing"
# note that i've chosen to use Hz instead of kHz

Plots.default(overwrite_figure = false) # causes display() to make a new window

function c_signal(x) # x is a time in seconds
	k_a = 1 # Hz
	k_b = 2 # Hz
	Φ = 3 * pi / 4
	return sin(2 * pi * k_a * x) + 1//2 * sin(2 * pi * k_b * x + Φ)
end

# signal as a collection of points
fᵢ = 500 # signal points frequency, Hz
Nᵢ = 1500 # number of samples for the signal data
tᵢ = 1 / fᵢ # signal points period, seconds
signal_points = [n * tᵢ for n in 0:(Nᵢ - 1)]
signal_data = [c_signal(x) for x in signal_points]

# note that this plots unitless-index vs value
# unitless index | 0:(Nᵢ - 1)
signal_plot = plot(signal_points, signal_data, title="signal", label="Nᵢ = $Nᵢ, kᵢ = $fᵢ")
#display(signal_plot)

fₛ = 8 # sampling frequency. samples per meter
Nₛ = 24 # number of samples
tₛ = 1 / fₛ # sampling period, seconds
sample_points = [n * tₛ for n in 0:(Nₛ - 1)]
sample_data = [c_signal(x) for x in sample_points]
# unitless index | 0:(Nₛ - 1)
sample_plot = scatter(sample_points, sample_data, title="sample", label="Nₛ = $Nₛ, kₛ = $fₛ")

# showing small table like lyons does
for n in 1:Nₛ
	println("	x($(n)) = $(sample_data[n]),")
end

# ℱ.signal
ℱ_signal_data = fft(signal_data)
# unitless index | 0:(Nᵢ - 1)
freq_signal_points = [(fᵢ / Nᵢ) * n for n in 0:(Nᵢ - 1)]
ℱ_signal_plot = plot(freq_signal_points, real.(ℱ_signal_data), title="ℱ.signal", label="real")

# ℱ.sample
ℱ_sample_data = fft(sample_data)
freq_sample_points = [(fₛ / Nₛ) * n for n in 0:(Nₛ - 1)]

# showing small table like lyons does
for n in 1:Nₛ
	println("	X($(n - 1)) = $(ℱ_sample_data[n]),")
end

ℱ_sample_plot = scatter(
	freq_sample_points,
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