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

I = ComplexF64.(make_field("in/$illumination_file_name.png"))
#save("out/$file_name--I.png", cap(abs.(I)))

Gˈ[0] = Complex{Float64}.(make_field("in/$(file_name).png")) # read in pattern. cast to Complex64 type.
image_size = size(Gˈ[0])[1]
random_phases = (rand(Float64, size(Gˈ[0])) .- 0.5) .* 2π
# ^ generate random phases. start with r ∈ (0 to 1]. then rescale to fit in (-π to π].
gˈ[0] = F⁻¹(Gˈ[0])# create relief guess by reverse-propagating pattern.
gˈ[0] = set_phase(gˈ[0], random_phases)
g[0] = set_modulus(gˈ[0], 1) # force-satisfy constraints for first object guess.

Hˈ[0] = copy(Gˈ[0])
hˈ[0] = copy(gˈ[0])
h[0] = copy(g[0])

β = 0.5
function fienup_fix(zˈ, z)
	if abs(zˈ) > 0
		return zˈ
	else
		return z - (β * zˈ)
	end
end

function scale_fix(zˈ, z)
	# abs(z) is 1, assumed, bc it satisfies object constraints.
	if abs(zˈ) < 1
		return zˈ + β * z
	elseif abs(zˈ) > 1
		return zˈ - β * z
	else
		return zˈ
	end
end

for i in 1:total_cycles
	#g[i] = broadcast(fienup_fix, gˈ[i - 1], g[i - 1])
	g[i] = scale_fix.(gˈ[i - 1], g[i - 1])
	#g[i] = set_modulus(gˈ[i - 1], gˈ[0])
	#g[i] = set_modulus(gˈ[i - 1], 1)
	G[i] = F(g[i])
	Gˈ[i] = set_modulus(G[i], Gˈ[0])
	gˈ[i] = F⁻¹(Gˈ[i])
	
	h[i] = set_modulus(hˈ[i - 1], hˈ[0])
	H[i] = F(h[i])
	Hˈ[i] = set_modulus(H[i], Hˈ[0]) # set the modulus of the current guess to be that of the target
	hˈ[i] = F⁻¹(Hˈ[i])

	#G_errors[i] = f_mse(abs.(G[i]), abs.(Gˈ[0]))
	G_errors[i] = f_mse(abs.(G[i]), abs.(Gˈ[i]))
	#G_errors[i] = norm(abs.(G[i]), abs.(Gˈ[0]))
	println("G_errors[$i] = ", □(G_errors[i]))
	#H_errors[i] = f_mse(abs.(H[i]), abs.(Hˈ[0]))
	H_errors[i] = f_mse(H[i], Hˈ[i])
	#H_errors[i] = norm(abs.(H[i]), abs.(Hˈ[0]))
	#println("H_errors[$i] = ", □(H_errors[i]))

	#png(heatmap(abs.(Gˈ[i]), title="|G′|"), "out/$(file_name)--moduli-G'$i.png")
	png(heatmap(abs.(half_roll(gˈ[i])), title="|g′|"), "out/$(file_name)--moduli-gg'$(i).png")
	#png(heatmap(angle.(Gˈ[i]), title="Φ(G′)"), "out/$(file_name)--phases-G'$i.png")
	png(heatmap(angle.(half_roll(gˈ[i])), title="Φ(g′)"), "out/$(file_name)--phases-gg'$(i).png")

	#png(heatmap(abs.(G[i]), title="|G|"), "out/$(file_name)--moduli-G$i.png")
	#png(heatmap(abs.(half_roll(g[i])), title="|g|"), "out/$(file_name)--moduli-gg$(i).png")
	#png(heatmap(angle.(half_roll(g[i])), title="Φ(g)"), "out/$(file_name)--phases-gg$(i).png")

	#png(histogram(vec(abs.(G[i])), bins=20, title="|G| distribution", yaxis = :log), "out/$(file_name)--moduli-G$i-distribution.png")
	#png(histogram(vec(angle.(g[i])), bins=20, title="Φ(g) distribution",), "out/$(file_name)--phases-gg$i-distribution.png")

	#save("out/$(file_name)--phases-gg$(i)-matrix.png", cap(angle.(half_roll(g[i]))))
	println("Cycle $i finished.")
end

gˈ★ = gˈ[total_cycles] # \bigstar
g★ = g[total_cycles] # \bigstar
save("out/$file_name--moduli-gg'-final.png", half_roll(cap(abs.(gˈ★))))
save("out/$file_name--moduli-G-final.png", cap(abs.(F(g★))))
save("out/$file_name--phases-gg-final.png", half_roll(cap(angle.(g★))))

garnet_color = RGB(((115, 0, 17) ./ 255)...)
orangish_color = RGB(((245, 102, 0) ./ 255)...)

#errors_max = max(maximum(G_errors), maximum(H_errors))
#G_errors_step = (errors_max - minimum(G_errors)) / 10#length(G_errors)
G_errors_step = (maximum(G_errors) - minimum(G_errors)) / length(G_errors)
errors_plot = plot(G_errors,
				   label = "G",
				   linecolor = false,
				   fill = (0, garnet_color))
#errors_plot = plot!(errors_plot, H_errors,
#				   label = "H",
#				   linecolor = false, # turns off line color
#				   fill = (0, orangish_color))
errors_plot = plot!(errors_plot,
					title = "Errors from Target Pattern",
				 	aspect_ratio = :auto, # VERY IMPORTANT
				 	xlabel = "Cycle",
				 	ylabel = "Mean Squared Error (MSE)",
				 	xlims = (1, length(G_errors)),
				 	ylims = (minimum(G_errors), maximum(G_errors)),
				 	yticks = ■(3).(minimum(G_errors):G_errors_step:maximum(G_errors)),) 
png(errors_plot, "out/$file_name--errors.png")

# filter off the high-frequency values. IMPORTANT: do not roll gˈ★[i] before doing this!
#gˈ_final_filtered = gˈ★[1:lffs, 1:lffs]
#Gˈ_final_filtered = F(vcat(hcat(gˈ_final_filtered, zeros(lffs, image_size - lffs)), zeros(image_size - lffs, image_size)))
#save("out/$file_name--G'-final--low-freq-$lffs.png", cap(abs.(Gˈ_final_filtered)))
#save("out/$file_name--g'-final--low-freq-$lffs.png", cap(abs.(gˈ_final_filtered)))

#png(heatmap(I, title="I"), "out/$(file_name)--I.png")
q = half_roll(g★) # basically q at least
#g_illuminated = set_modulus(q, abs.(I) .* abs.(q))
g_illuminated = set_modulus(q, abs.(I))
#g_illuminated = I .* (Φ	-> exp(1im * Φ)).(angle.(q))
save("out/$file_name--gg-illuminated.png", cap(abs.(g_illuminated)))
save("out/$file_name--G-illuminated.png", cap(abs.(F(g_illuminated))))

end # end module