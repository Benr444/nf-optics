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

G[0] = OpticField("in/three.png")
g[0] = F⁻¹(G[0])
png(g[0], "out/three-g0.png")
total_cycles = 6
for i in 1:total_cycles
	gˈ[i] = set_phase(g[i - 1], 0)
	G[i] = F(gˈ[i])
	Gˈ[i] = set_modulus(G[i], G[0].values)
	g[i] = F⁻¹(Gˈ[i])
	println("Cycle $i finished.")
	#r_Gˈ = (abs.(G[i].values))
	png(heatmap(
		heatmap(abs.(gˈ[i].values), title="|gˈ|"),
		heatmap(abs.(G[i].values), title="|G|"),
		heatmap(abs.(Gˈ[i].values), title="|Gˈ|"),
		heatmap(abs.(g[i].values), title="|g|"),
		layout=4, plot_title="absolute values"), "out/three-moduli$i.png")
	png(heatmap(
		heatmap(angle.(gˈ[i].values), title="Φ(gˈ)"),
		heatmap(angle.(G[i].values), title="Φ(G)"),
		heatmap(angle.(Gˈ[i].values), title="Φ(Gˈ)"),
		heatmap(angle.(g[i].values), title="Φ(g)"),
		layout=4, plot_title="phase"), "out/three-phases$i.png")
end

end # end module