# usage
# from cmd: 
# -> "julia" 
# -> "include("DesignKinoform.jl")"
# -> re-run the previous after every edit.
# for best performance, do not close julia\

# note that the 'angle' function returns an value in the range (-π to π).

module DesignKinoform

using Revise
using Images, FileIO, FFTW, Random
using Plots; pyplot()
Revise.includet("ROptics.jl")
using Main.ROptics

file_name = "uofsc-200px"

Plots.default(aspect_ratio = :auto, size = (800, 800))

my_λ = 1.0 # units:
my_scale = 1.0 # units:
low_freq_filter_size = lffs = 20
total_cycles = 20

g, gˈ, G, Gˈ = Dict(), Dict(), Dict(), Dict() # prime marker is \verts
# g = relief guess.
# gˈ = phase-zero'd relief guess (same moduli as g, since that is unaffected).
# G = far-field pattern.
# Gˈ = moduli-corrected pattern (same moduli as given initial pattern).

h, hˈ, H, Hˈ = Dict(), Dict(), Dict(), Dict()
# binary-level data

G_errors = zeros(total_cycles)
# ^ stores diffs from the target G[0]. there is 1 such diff for each cycle (0 is not graphed).
H_errors = zeros(total_cycles)

G[0] = Complex{Float64}.(make_field("in/$(file_name).png")) # read in pattern. cast to Complex64 type.
image_size = size(G[0])[1]
#G[0] = set_phase(G[0], 0)
g[0] = F⁻¹(G[0]) # create relief guess by reverse-propagating pattern
random_phases = (rand(Float64, size(g[0])) .- 0.5) .* 2π
# ^ generate random phases. start with r ∈ (0 to 1]. then rescale to fit in (-π to π].
g[0] = set_phase(g[0], random_phases)
#png(heatmap(abs.(G[0]), title="|G|"), "out/$(file_name)--moduli-G0.png")
#png(heatmap(abs.(half_roll(g[0])), title="|g|"), "out/$(file_name)--moduli-gg0.png")
#png(heatmap(angle.(G[0]), title="Φ(G)"), "out/$(file_name)--phases-G0.png")
#png(heatmap(angle.(half_roll(g[0])), title="Φ(g)"), "out/$(file_name)--phases-gg0.png")

H[0] = copy(G[0])
h[0] = copy(g[0])

for i in 1:total_cycles
	#gˈ[i] = set_phase(g[i - 1], 0)
	gˈ[i] = set_modulus(g[i - 1], g[0])
	G[i] = F(gˈ[i])
	Gˈ[i] = set_modulus(G[i], G[0]) # set the modulus of the current guess to be that of the target
	g[i] = F⁻¹(Gˈ[i])

	# cast so hˈ[i] only has two values
	# in the end, it will be rescaled to maximally fit in -π to \pi
	binary_phases = (c -> angle(c) < 0 ? -π : π).(h[i - 1])
	hˈ[i] = set_phase(set_modulus(h[i - 1], h[0]), binary_phases)
	H[i] = F(hˈ[i])
	Hˈ[i] = set_modulus(H[i], H[0])
	h[i] = F⁻¹(Hˈ[i])
	
	G_errors[i] = mean_se(abs.(G[0]), abs.(G[i]))
	#println("G_errors[$i] = ", □(G_errors[i]))
	H_errors[i] = mean_se(abs.(H[0]), abs.(H[i]))
	
	#png(heatmap(abs.(G[i]), title="|G|"), "out/$(file_name)--moduli-G$i.png")
	#png(heatmap(abs.(half_roll(g[i])), title="|g|"), "out/$(file_name)--moduli-gg$(i).png")
	#png(heatmap(angle.(G[i]), title="Φ(G)"), "out/$(file_name)--phases-G$i.png")
	#png(heatmap(angle.(half_roll(g[i])), title="Φ(g)"), "out/$(file_name)--phases-gg$(i).png")

	#png(histogram(vec(abs.(G[i])), bins=20, title="|G| distribution", yaxis = :log), "out/$(file_name)--moduli-G$i-distribution.png")
	png(histogram(vec(angle.(g[i])), bins=20, title="Φ(g) distribution", yaxis = :log), "out/$(file_name)--phases-gg$i-distribution.png")

	#save("out/$(file_name)--phases-gg$(i)-matrix.png", cap(angle.(half_roll(g[i]))))
	println("Cycle $i finished.")
end

save("out/$(file_name)--moduli-G-final-matrix.png", cap(abs.(G[total_cycles])))

garnet_color = RGB(((115, 0, 17) ./ 255)...)
orangish_color = RGB(((245, 102, 0) ./ 255)...)

@show H_errors
errors_max = max(maximum(G_errors), maximum(H_errors))

G_errors_step = (errors_max - minimum(G_errors)) / 10#length(G_errors)
errors_plot = plot(H_errors,
				   label = "Quantized object",
				   linecolor = false, # turns off line color
				   fill = (0, orangish_color))
errors_plot = plot!(errors_plot, G_errors,
				    label = "Continuous object",
				    linecolor = false,
				    fill = (0, garnet_color))
errors_plot = plot!(errors_plot,
					title = "Errors from Target Pattern",
				 	aspect_ratio = :auto, # VERY IMPORTANT
				 	xlabel = "Cycle",
				 	ylabel = "MSE",
				 	xlims = (1, length(G_errors)),
				 	ylims = (minimum(G_errors), maximum(G_errors)),
				 	yticks = ■(3).(minimum(G_errors):G_errors_step:maximum(G_errors)),) 

png(errors_plot, "out/$file_name--errors.png")

# filter off the high-frequency values. IMPORTANT: do not roll g[i] before doing this!
g_final_filtered = g[total_cycles][1:lffs, 1:lffs]
G_final_filtered = F(vcat(hcat(g_final_filtered, zeros(lffs, image_size - lffs)), zeros(image_size - lffs, image_size)))

save("out/$file_name--G-final--low-freq-$lffs.png", cap(abs.(G_final_filtered)))
save("out/$file_name--gg-final.png", half_roll(cap(abs.(g[total_cycles]))))
save("out/$file_name--gg-final--low-freq-$lffs.png", cap(abs.(g_final_filtered)))

end # end module