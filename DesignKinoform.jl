# usage
# from cmd: 
# -> "julia" 
# -> "include("DesignKinoform.jl")"
# -> re-run the previous after every edit.
# for best performance, do not close julia\

# note that the 'angle' function returns an value in the range (-π to π).

module DesignKinoform

using Revise
using Images, FileIO, Plots, FFTW
Revise.includet("ROptics.jl")
using Main.ROptics

file_name = "three-200px"

Plots.default(aspect_ratio = 1, size = (800, 800))

my_λ = 1.0 # units:
my_scale = 1.0 # units:
g, gˈ, G, Gˈ = Dict(), Dict(), Dict(), Dict() # prime marker is \verts
# g = relief guess.
# gˈ = phase-zero'd relief guess (same moduli as g, since that is unaffected).
# G = far-field pattern.
# Gˈ = moduli-corrected pattern (same moduli as given initial pattern).

G[0] = Complex{Float64}.(make_field("in/$(file_name).png")) # read in pattern. cast to Complex64 type.
g[0] = F⁻¹(G[0]) # create relief guess by reverse-propagating pattern
png(heatmap(abs.(G[0]), title="|G|"), "out/$(file_name)--moduli-G0.png")
png(heatmap(abs.(half_roll(g[0])), title="|g|"), "out/$(file_name)--moduli-gg0.png")
total_cycles = 6
for i in 1:total_cycles
	gˈ[i] = set_phase(g[i - 1], 0)
	G[i] = F(gˈ[i])
	Gˈ[i] = set_modulus(G[i], G[0]) # set the modulus of the current guess to be that of the target
	g[i] = F⁻¹(Gˈ[i])
	println("Cycle $i finished.")
	png(heatmap(abs.(G[i]), title="|G|"), "out/$(file_name)--moduli-G$i.png")
	png(heatmap(abs.(half_roll(g[i])), title="|g|"), "out/$(file_name)--moduli-gg$(i).png")
	#png(heatmap(angle.(G[i]), title="Φ(G)"), "out/$(file_name)--phases-G$i.png")
	#png(heatmap(angle.(half_roll(g[i])), title="Φ(g)"), "out/$(file_name)--phases-gg$(i).png")

	Gi_moduli_list = vec(abs.(G[i]))
	#@show typeof(Gi_moduli_list)
	png(histogram(Gi_moduli_list, bins=20, title="|G| distribution", yaxis = :log), "out/$(file_name)--moduli-G$i-distribution.png")

	save("out/$(file_name)--moduli-gg$(i)-matrix.png", cap(abs.(half_roll(g[i]))))
end
	
super_matrix_row_1 = hcat(cap(abs.(half_roll(g[0]))), cap(abs.(half_roll(g[1]))))
super_matrix_row_2 = hcat(cap(abs.(half_roll(g[2]))), cap(abs.(half_roll(g[3]))))
super_matrix_row_3 = hcat(cap(abs.(half_roll(g[4]))), cap(abs.(half_roll(g[5]))))
super_matrix = vcat(super_matrix_row_1, super_matrix_row_2, super_matrix_row_3)
save("out/$(file_name)--super-matrix.png", super_matrix)

end # end module