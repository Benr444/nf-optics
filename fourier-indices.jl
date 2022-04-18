using Revise
using Images, FileIO, Plots, FFTW

# https://jackschaedler.github.io/circles-sines-signals/dft_frequency.html

phasor(f, t) =	exp(2 * pi * 1im * f * t)

samples = [0.0 + 0.0im,
		   2.3 + 3.3im,
		   4.4 + 12.3im,
		   -3.1 + -32.0im,
		   430.0 + -0.01im]

function get_Amplitude(f, data)
	length(data)
	duration = 1/f

	for t_index in 1:10
		data[] * phasor(-f, t)
	end
end