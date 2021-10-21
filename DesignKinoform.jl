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

my_λ = 1.0 # units:
my_scale = 1.0 # units:
g, gˈ, G, Gˈ = Dict(), Dict(), Dict(), Dict() # prime marker is \verts

G[0] = make_field("in/three.png")
g[0] = F⁻¹(G[0])
png(g[0], "out/three-g0.png")
total_cycles = 6
for i in 1:total_cycles
	gˈ[i] = set_phase(g[i - 1], 0)
	G[i] = F(gˈ[i])
	Gˈ[i] = set_modulus(G[i], G[0])
	g[i] = F⁻¹(Gˈ[i])
	println("Cycle $i finished.")
	#r_Gˈ = (abs.(G[i]))
	png(heatmap(
		heatmap(abs.(gˈ[i]), title="|gˈ|"),
		heatmap(abs.(G[i]), title="|G|"),
		heatmap(abs.(Gˈ[i]), title="|Gˈ|"),
		heatmap(abs.(g[i]), title="|g|"),
		layout=4, plot_title="absolute values"), "out/three-moduli$i.png")
	png(heatmap(
		heatmap(angle.(gˈ[i]), title="Φ(gˈ)"),
		heatmap(angle.(G[i]), title="Φ(G)"),
		heatmap(angle.(Gˈ[i]), title="Φ(Gˈ)"),
		heatmap(angle.(g[i]), title="Φ(g)"),
		layout=4, plot_title="phase"), "out/three-phases$i.png")
end

end # end module