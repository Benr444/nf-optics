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
total_cycles = 50

g, gˈ, G, Gˈ = Dict(), Dict(), Dict(), Dict() # prime marker is \verts
# g = relief guess.
# gˈ = phase-zero'd relief guess (same moduli as g, since that is unaffected).
# G = far-field pattern.
# Gˈ = moduli-corrected pattern (same moduli as given initial pattern).

G_errors = zeros(total_cycles)
# ^ stores diffs from the target G[0]. there is 1 such diff for each cycle (0 is not graphed).

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
for i in 1:total_cycles
	#gˈ[i] = set_phase(g[i - 1], 0)
	gˈ[i] = set_modulus(g[i - 1], g[0])
	G[i] = F(gˈ[i])
	Gˈ[i] = set_modulus(G[i], G[0]) # set the modulus of the current guess to be that of the target
	g[i] = F⁻¹(Gˈ[i])

	#@show gˈ_modulus_mean = mean(abs.(gˈ[i]))
	#@show gˈmodulus_diff = round(norm(abs.(gˈ[i]), 1.0), sigdigits = 3)
	#@show gˈ_modulus_is_1 = all(abs.(gˈ[i]) .== 1.0)

	#@show sizes_match = (size(g[i]) == size(G[i]))

	G_raw = F(g[i - 1]) # uncorrected propagation.
	#@show correction_diff = round(norm(abs.(G_raw), abs.(G[i])), sigdigits = 3)

	#@show insanity_modulus_diff = round(norm(abs.(G_raw), abs.(G[0])), sigdigits = 3)
	#@show insanity_phase_diff = round(norm(angle.(G_raw), angle.(G[0])), sigdigits = 3)
	
	#G_errors[i] = norm(cap(abs.(G[0])), cap(abs.(G[i])))
	G_errors[i] = norm(abs.(G[0]), abs.(G[i]))
	#@show round(G_errors[i], sigdigits = 3)
	println("G_errors[$i] = ", □(G_errors[i]))
	#@show norm(abs.(G[0])) |> □
	#@show norm_G_i = norm(abs.(G[i]))
	#@show norm_gˈ_i = norm(abs.(gˈ[i]))
	#@show plancherel_diff = norm_gˈ_i - norm_G_i
	#@show size(gˈ[i]) == size(G[i])
	#@show Gˈ_error = round(norm(abs.(G[0]), abs.(Gˈ[i])), sigdigits = 3)
	#@show g_error = round(norm(abs.(g[0]), abs.(g[i])), sigdigits = 3)
	#@show gˈ_error = round(norm(abs.(g[0]), abs.(gˈ[i])), sigdigits = 3)

	#@show g_change = round(norm(abs.(g[i - 1]), abs.(g[i])), sigdigits = 3)

	png(heatmap(abs.(G[i]), title="|G|"), "out/$(file_name)--moduli-G$i.png")
	png(heatmap(abs.(half_roll(g[i])), title="|g|"), "out/$(file_name)--moduli-gg$(i).png")
	png(heatmap(angle.(G[i]), title="Φ(G)"), "out/$(file_name)--phases-G$i.png")
	png(heatmap(angle.(half_roll(g[i])), title="Φ(g)"), "out/$(file_name)--phases-gg$(i).png")

	#Gi_moduli_list = vec(abs.(G[i]))
	#png(histogram(Gi_moduli_list, bins=20, title="|G| distribution", yaxis = :log), "out/$(file_name)--moduli-G$i-distribution.png")

	#save("out/$(file_name)--phases-gg$(i)-matrix.png", cap(angle.(half_roll(g[i]))))
	println("Cycle $i finished.")
end


save("out/$(file_name)--moduli-G-final-matrix.png", cap(abs.(G[total_cycles])))

G_errors_step = (maximum(G_errors) - minimum(G_errors)) / 10#length(G_errors)
G_errors_plot = plot(G_errors,
					 title = "Errors of G",
					 label = "",
					 aspect_ratio = :auto, # VERY IMPORTANT
					 xlabel = "Cycle",
					 ylabel = "Sum Square Differences of Moduli of G relative to G[0]",
					 xlims = (1, length(G_errors)),
					 #series_annotations = string.(round.(G_errors, sigdigits = 3)),
					 ylims = (minimum(G_errors), maximum(G_errors)),
					 yticks = (x -> round(x, sigdigits = 5)).(minimum(G_errors):G_errors_step:maximum(G_errors)))

png(G_errors_plot, "out/$file_name--G-errors.png")

# filter off the high-frequency values. IMPORTANT: do not roll g[i] before doing this!
g_final_filtered = g[total_cycles][1:lffs, 1:lffs]
G_final_filtered = F(vcat(hcat(g_final_filtered, zeros(lffs, image_size - lffs)), zeros(image_size - lffs, image_size)))

@show size(cap(abs.(G[total_cycles])))
@show size(cap(abs.(G_final_filtered)))

save("out/$file_name--G-final--low-freq-$lffs.png", cap(abs.(G_final_filtered)))
save("out/$file_name--gg-final.png", half_roll(cap(abs.(g[total_cycles]))))
save("out/$file_name--gg-final--low-freq-$lffs.png", cap(abs.(g_final_filtered)))

end # end module