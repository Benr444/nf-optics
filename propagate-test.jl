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

G[0] = Complex{Float64}.(make_field("in/$(file_name).png")) # read in pattern
g[0] = F⁻¹(G[0]) # create relief guess by reverse-propagating pattern
png(heatmap(abs.(G[0]), title="|G|"), "out/$(file_name)-moduli-G0.png")
#png(heatmap(abs.(half_roll(g[0])), title="|g|"), "out/$(file_name)-moduli-lilg0.png")

G0_recreation = F(F⁻¹(G[0]))
png(heatmap(abs.(G0_recreation), title="|G|"), "out/$(file_name)-moduli-G0-recreation.png")

ΔG0 = broadcast(-, G[0], G0_recreation)

png(heatmap(abs.(ΔG0), title="|G|"), "out/$(file_name)-moduli-G0-diff.png")
png(heatmap(angle.(ΔG0), title="|G|"), "out/$(file_name)-phases-G0-diff.png")

@show matching_moduli = all(broadcast(==, abs.(ΔG0), 0))
@show mean_moduli_diff = sum(abs.(ΔG0)) / length(ΔG0)
@show matching_phases = all(broadcast(==, angle.(ΔG0), 0))
@show mean_phase_diff = sum(abs.(angle.(ΔG0))) / length(ΔG0)

#png(heatmap(abs.(G[i]), title="|G|"), "out/$(file_name)-moduli-G$i.png")
#png(heatmap(abs.(half_roll(g[i])), title="|g|"), "out/$(file_name)-moduli-lilg$(i).png")
#png(heatmap(angle.(G[i]), title="Φ(G)"), "out/$(file_name)-phases-G$i.png")
#png(heatmap(angle.(half_roll(g[i])), title="Φ(g)"), "out/$(file_name)-phases-lilg$(i).png")

end # end module