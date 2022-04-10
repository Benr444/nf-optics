# usage
# from cmd: 
# -> "julia" 
# -> "include("DesignKinoform.jl")"
# -> re-run the previous after every edit.
# for best performance, do not close julia

module DesignKinoform

using Revise
using Images, FileIO, Plots, FFTW
Revise.includet("ROptics.jl")
using Main.ROptics

file_name = "three-200px"


my_λ = 1.0 # units:
my_scale = 1.0 # units:
g, gˈ, G, Gˈ = Dict(), Dict(), Dict(), Dict() # prime marker is \verts
# g = relief guess.
# gˈ = phase-zero'd relief guess (same moduli as g, since that is unaffected).
# G = far-field pattern.
# Gˈ = moduli-corrected pattern (same moduli as given initial pattern).

G[0] = make_field("in/$(file_name).png") # read in pattern
g[0] = F⁻¹(G[0]) # create relief guess by reverse-propagating pattern
png(heatmap(abs.(G[0]), title="|G|"), "out/$(file_name)-moduli-G0.png")
png(heatmap(abs.(half_roll(g[0])), title="|g|"), "out/$(file_name)-moduli-lilg0.png")
total_cycles = 6
for i in 1:total_cycles
	gˈ[i] = set_phase(g[i - 1], 0)
	G[i] = F(gˈ[i])
	Gˈ[i] = set_modulus(G[i], G[0]) # set the modulus of the current guess to be that of the target
	g[i] = F⁻¹(Gˈ[i])
	println("Cycle $i finished.")
	#r_Gˈ = (abs.(G[i]))
	#png(heatmap(
	#	#heatmap(abs.(gˈ[i]), title="|gˈ|"), # no point in showing this, will be same as g.
	#	heatmap(abs.(G[i]), title="|G|"),
	#	#heatmap(abs.(Gˈ[i]), title="|Gˈ|"), # no point in showing this, will be same as G[0].
	#	heatmap(abs.(g[i]), title="|g|"),
	#	layout=2, plot_title="moduli"), "out/$(file_name)-moduli$i.png")

	png(heatmap(abs.(G[i]), title="|G|"), "out/$(file_name)-moduli-G$i.png")
	
	png(heatmap(abs.(half_roll(g[i])), title="|g|"), "out/$(file_name)-moduli-lilg$(i).png")
	
	png(heatmap(angle.(G[i]), title="Φ(G)"), "out/$(file_name)-phases-G$i.png")
	
	png(heatmap(angle.(half_roll(g[i])), title="Φ(g)"), "out/$(file_name)-phases-lilg$(i).png")
	
	#png(heatmap(
	#	#heatmap(angle.(gˈ[i]), title="Φ(gˈ)"), # no point in showing this, will be 0.
	#	heatmap(angle.(G[i]), title="Φ(G)"),
	#	#heatmap(angle.(Gˈ[i]), title="Φ(Gˈ)"), # no point in showing this, will be same as G.
	#	heatmap(angle.(g[i]), title="Φ(g)"),
	#	layout=2, plot_title="phase"), "out/$(file_name)-phases$i.png")
end

end # end module