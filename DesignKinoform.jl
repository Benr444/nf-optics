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
using Plots; gr()
Revise.includet("ROptics.jl")
using Main.ROptics

file_name = "uofsc-200px"
illumination_file_name = "illumination-200px"

Plots.default(aspect_ratio = :auto, size = (1200, 1200), fontfamily = font("Computer Modern"))

low_freq_filter_size = lffs = 150
total_cycles = 30

g, g′, G, G′ = Dict(), Dict(), Dict(), Dict() # prime marker is \prime
# g = relief guess.
# g′ = phase-zero'd relief guess (same moduli as g, since that is unaffected).
# G = far-field pattern.
# G′ = moduli-corrected pattern (same moduli as given initial pattern).

G′[0] = Complex{Float64}.(make_field("in/$(file_name).png")) # read in pattern. cast to Complex64 type.
pattern_size = size(G′[0])[1]

px_scale = 1 # units: micrometers. distance across one pixel edge.
# pos of cell 1 -> 0 (initial term)
# pos of cell 2 -> 1 * px_scale
# pos of cell i -> (i-1) * px_scale 
# freq of cell 1 -> 0 (constant term)
# freq of cell 2 -> 1 cycle per set, eg (length(G) * px_scale)^-1
# freq of cell i -> i-1 cycles per set
pattern_length = (pattern_size) * px_scale
pattern_base_tone = (pattern_size * px_scale)^-1
phys_pos(i) = (i - 1) * px_scale 
phys_pos(i, j) = sqrt(phys_pos(i)^2 + phys_pos(j)^2)
phys_freq(i) = (i - 1) * (pattern_base_tone)
phys_freq(i, j) = sqrt(phys_freq(i)^2 + phys_freq(j)^2)

h, h′, H, H′ = Dict(), Dict(), Dict(), Dict()
# binary-level data

G_errors = zeros(total_cycles)
# ^ stores diffs from the target G[0]. there is 1 such diff for each cycle (0 is not graphed).
g′_errors = zeros(total_cycles)
H_errors = zeros(total_cycles)

I = ComplexF64.(make_field("in/$illumination_file_name.png"))
#save("out/$file_name--I.png", cap(abs.(I)))


random_phases = (rand(Float64, size(G′[0])) .- 0.5) .* 2π
# ^ generate random phases. start with r ∈ (0 to 1]. then rescale to fit in (-π to π].
g′[0] = F⁻¹(G′[0])# create relief guess by reverse-propagating pattern.
g′[0] = set_phase(g′[0], random_phases)
g[0] = set_modulus(g′[0], 1) # force-satisfy constraints for first object guess.

H′[0] = copy(G′[0])
h′[0] = copy(g′[0])
h[0] = copy(g[0])

β = 0.9
function Δg_d(z′, z)
	if abs(z′) == 1.0
		return 0
	else
		return abs(1.0) * (z′ / abs(z′)) - z′
	end
end

for i in 1:total_cycles
	
	#g[i] = set_modulus(g′[i - 1], g′[0])
	#g[i] = set_modulus(g′[i - 1], mean(g′[i - 1]))
	g[i] = set_modulus(g′[i - 1], 1) # <----- this is really the one i swear. eq 24 in Fienup 1984 reduces to this if β = 1.
	#g[i] = set_modulus(g′[i - 1], β) .+ (1 - β) .* h′[i - 1]
	# this and the one below are algebraically equivalent (if i havent fucked with Δg_d since writing this comment).
	#g[i] = g′[i - 1] .+ (β .* Δg_d.(g′[i - 1], g[i - 1]))

	G[i] = F(g[i])
	G′[i] = set_modulus(G[i], G′[0])
	g′[i] = F⁻¹(G′[i])
	
	#h[i] = set_modulus(h′[i - 1], h′[0])
	h[i] = set_modulus(h′[i - 1], 1.0)
	H[i] = F(h[i])
	H′[i] = set_modulus(H[i], H′[0]) # set the modulus of the current guess to be that of the target
	h′[i] = F⁻¹(H′[i])

	G_errors[i] = f_mse(abs.(G[i]), abs.(G′[0]))
	#G_errors[i] = f_mse(abs.(G[i]), abs.(G′[i]))
	#g′_errors[i] = f_mse(abs.(g[i]), abs.(g′[i]))
	#G_errors[i] = norm(abs.(G[i]), abs.(G′[0]))
	println("G_errors[$i] = ", □(G_errors[i]))
	#println("g'_errors[$i] = ", □(g′_errors[i]))
	H_errors[i] = f_mse(abs.(H[i]), abs.(H′[0]))
	#H_errors[i] = f_mse(abs.(H[i]), abs.(H′[i]))
	#H_errors[i] = norm(abs.(H[i]), abs.(H′[0]))
	println("H_errors[$i] = ", □(H_errors[i]))

	#png(heatmap(abs.(G′[i]), title="|G′|"), "out/$(file_name)--moduli-G'$i.png")
	#png(heatmap(abs.(half_roll(g′[i])), title="|g′|"), "out/$(file_name)--moduli-gg'$(i).png")
	#png(heatmap(abs.(half_roll(g[i])), title="|g|"), "out/$(file_name)--moduli-gg$(i).png")
	#png(heatmap(angle.(G′[i]), title="Φ(G′)"), "out/$(file_name)--phases-G'$i.png")
	#png(heatmap(angle.(half_roll(g′[i])), title="Φ(g′)"), "out/$(file_name)--phases-gg'$(i).png")

	#png(heatmap(abs.(G[i]), title="|G|"), "out/$(file_name)--moduli-G$i.png")
	#png(heatmap(abs.(half_roll(g[i])), title="|g|"), "out/$(file_name)--moduli-gg$(i).png")
	#png(heatmap(angle.(half_roll(g[i])), title="Φ(g)"), "out/$(file_name)--phases-gg$(i).png")

	#png(histogram(vec(abs.(G[i])), bins=20, title="|G| distribution", yaxis = :log), "out/$(file_name)--moduli-G$i-distribution.png")
	#png(histogram(vec(angle.(g[i])), bins=20, title="Φ(g) distribution",), "out/$(file_name)--phases-gg$i-distribution.png")

	#save("out/$(file_name)--phases-gg$(i)-matrix.png", cap(angle.(half_roll(g[i]))))
	println("Cycle $i finished.")
end

g′★ = g′[total_cycles] # \bigstar
g★ = g[total_cycles] # \bigstar
G′★ = G′[total_cycles] # \bigstar
G★ = G[total_cycles] # \bigstar

save("out/$file_name--moduli-gg-final.png", half_roll(cap(abs.(g★))))
save("out/$file_name--moduli-gg'-final.png", half_roll(cap(abs.(g′★))))
save("out/$file_name--moduli-G-final.png", cap(abs.(F(g★))))
save("out/$file_name--phases-gg-final.png", half_roll(cap(angle.(g★))))

#n = size(g′★)[1]
#@show plancheral_G′ = sum(abs.(G′★).^2)
#@show plancheral_g′ = sum(abs.(g′★).^2)
#@show (n^2) * plancheral_g′ - plancheral_G′

garnet_color = RGB(((115, 0, 17) ./ 255)...)
orangish_color = RGB(((245, 102, 0) ./ 255)...)

#errors_max = max(maximum(G_errors), maximum(H_errors))
#G_errors_step = (errors_max - minimum(G_errors)) / 10#length(G_errors)
G_errors_step = (maximum(G_errors) - minimum(G_errors)) / 10#length(G_errors)
errors_plot = plot(G_errors,
				   label = "G",
				   linecolor = false,
				   fill = (0, garnet_color))
errors_plot = plot!(errors_plot, H_errors,
				   label = "H",
				   linecolor = false, # turns off line color
				   fill = (0, orangish_color))
errors_plot = plot!(errors_plot,
					title = "Errors from Target Pattern",
				 	aspect_ratio = :auto, # VERY IMPORTANT
				 	xlabel = "Cycle",
				 	ylabel = "Mean Squared Error (MSE)",
				 	xlims = (1, length(G_errors)),
				 	ylims = (minimum(G_errors), maximum(G_errors)),
				 	yticks = ■(3).(minimum(G_errors):G_errors_step:maximum(G_errors)),) 
png(errors_plot, "out/$file_name--errors.png")

# filter off the high-frequency values. IMPORTANT: do not roll g′★[i] before doing this!
g★_lows = g★[1:lffs, 1:lffs]
g★_filtered = vcat(hcat(g★_lows, zeros(lffs, pattern_size - lffs)), zeros(pattern_size - lffs, pattern_size))
G★_filtered = F(g★_filtered)
save("out/$file_name--G-star--low-freq-$lffs.png", cap(abs.(G★_filtered)))
png(heatmap(abs.(g★_lows)), "out/$file_name--gg-star--low-freq-$lffs.png")

#png(heatmap(I, title="I"), "out/$(file_name)--I.png")
q = half_roll(g★) # basically q at least
#g_illuminated = set_modulus(q, abs.(I) .* abs.(q))
g_illuminated = set_modulus(q, abs.(I))
#g_illuminated = I .* (Φ	-> exp(1im * Φ)).(angle.(q))
save("out/$file_name--gg-illuminated.png", cap(abs.(g_illuminated)))
save("out/$file_name--G-illuminated.png", cap(abs.(F(g_illuminated))))

illumination_error = mean_se(abs.(F(g_illuminated)), abs.(G′[total_cycles]))
println("illumination_error = ", illumination_error)

end # end module