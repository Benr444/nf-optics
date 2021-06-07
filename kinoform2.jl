### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 02ea27a9-aef2-4f9f-91fe-898ab0236900
using Images, FileIO, Plots, FFTW

# ╔═╡ 5aeab850-c7a7-11eb-2b60-b50db938f0dd
md"# kinoform"

# ╔═╡ 84e0c901-0974-4361-90b9-234f924fa254
import Base.real

# ╔═╡ 55f41df0-8996-477b-a912-b4d17b076e81
struct OpticParams
	λ::Float64
end

# ╔═╡ 6cb0d982-f1dd-4897-ad6b-840db070878d
# OpticParams() = OpticParams(1.0)

# ╔═╡ 45be7f7e-f579-4e56-b5bd-21c07ea52426
struct OpticField # the minimum data for an optic-design
	values::Matrix{Complex{Float64}} # 2D Array of Complex numbers
	params::OpticParams
	#scale::Float64 # physical-distance between matrix slots
	#λ::Float64 # physical wavelength
end

# ╔═╡ 5dc92631-902f-48b9-a8d9-ce62093d6c27
make_OpticField(depths::Matrix{Complex{Float64}}; params::OpticParams = OpticParams(1.0)) = OpticField(depths, params)

# ╔═╡ 3d2f604a-6ff6-4b02-8976-eec51986e687
shape(x::OpticField) = shape(x.values)

# ╔═╡ e551e745-209d-47ff-bfab-0d6a05611998
set_phase(z, Φ) = complex(abs(z) * cos(Φ), abs(z) * sin(Φ))

# ╔═╡ ed16009e-8b5c-4181-a344-65ed03993fe1
set_modulus(z, l) = complex(abs(l) * cos(angle(z)), abs(l) * sin(angle(z)))

# ╔═╡ 6e1ea2d0-6dd2-47ce-893d-1da75f14abc1
function grayModulus(c::AbstractRGBA)
	colors = (red(c), green(c), blue(c))
	(sqrt ∘ sum)((colors.^2))
end

# ╔═╡ 6c50252b-8f91-439c-a808-e92fa5ff5d14
function make_OpticField(filename::String; params::OpticParams = OpticParams(1.0))
	image = load(filename)
	gray_array = grayModulus.(image)
	OpticField(gray_array, params)
end

# ╔═╡ 0e6c3c01-bb0f-4125-ba8e-30a211f81fbc
F(x::OpticField) = make_OpticField(fft(x.values); x.params) # forward-propagate

# ╔═╡ c3035153-03e4-4459-bd8e-3575e61ee2b2
F⁻¹(x::OpticField) = make_OpticField(ifft(x.values), x.params) # reverse-propagate

# ╔═╡ 897edd6f-a9f8-4da4-ac34-031445495752
set_phase(x::OpticField, Φs) = make_OpticField(broadcast(set_phase, x.values, Φs), x.params)

# ╔═╡ 74d4b1a2-91d9-45d9-b605-3113b9268794
set_modulus(x::OpticField, Φs) = make_OpticField(broadcast(set_modulus, x.values, Φs), x.params)

# ╔═╡ 6d2c0ec7-4b4a-409a-98e5-dc16c362558a
function fitMyRange(xx)
	xx_min = minimum(xx)
	range = maximum(xx) - xx_min
	(xx .- xx_min)./range # -1 -> 1
end

# ╔═╡ 6e1c99b0-8d89-48d7-9226-9c98690df00d
Plots.png(x::OpticField, filepath) = png(plot(Gray24.(fitMyRange(real.(x.values)))), filepath)

# ╔═╡ 8d9a068b-42aa-441c-86f8-7beab274f27c
my_λ = 1.0 # units:

# ╔═╡ 49546fa4-274d-4a62-8138-05046bfc7450
my_scale = 1.0 # units:

# ╔═╡ d86e29cb-b67d-41b5-bb4a-9199a1e8b3db
g, gˈ, G, Gˈ = Dict(), Dict(), Dict(), Dict() # prime marker is \verts

# ╔═╡ e45a2654-1a19-49c4-8a44-7b859dc0f5d1
G[0] = make_OpticField("in/three.png")

# ╔═╡ 0923a7e7-4838-4d6a-babf-6551c1b0fa5d
g[0] = F⁻¹(G[0])

# ╔═╡ 9d3afdde-7903-4539-ac3e-b56c38c4c834
png(g[0], "out/three-g0.png")

# ╔═╡ f539ff24-56cf-473d-a29f-66fafa579d50
total_cycles = 3

# ╔═╡ abbe0b94-f93c-4614-a534-020fc1526301
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

# ╔═╡ Cell order:
# ╠═5aeab850-c7a7-11eb-2b60-b50db938f0dd
# ╠═02ea27a9-aef2-4f9f-91fe-898ab0236900
# ╠═84e0c901-0974-4361-90b9-234f924fa254
# ╠═55f41df0-8996-477b-a912-b4d17b076e81
# ╠═6cb0d982-f1dd-4897-ad6b-840db070878d
# ╠═45be7f7e-f579-4e56-b5bd-21c07ea52426
# ╠═5dc92631-902f-48b9-a8d9-ce62093d6c27
# ╠═6c50252b-8f91-439c-a808-e92fa5ff5d14
# ╠═3d2f604a-6ff6-4b02-8976-eec51986e687
# ╠═0e6c3c01-bb0f-4125-ba8e-30a211f81fbc
# ╠═c3035153-03e4-4459-bd8e-3575e61ee2b2
# ╠═6e1c99b0-8d89-48d7-9226-9c98690df00d
# ╠═e551e745-209d-47ff-bfab-0d6a05611998
# ╠═897edd6f-a9f8-4da4-ac34-031445495752
# ╠═ed16009e-8b5c-4181-a344-65ed03993fe1
# ╠═74d4b1a2-91d9-45d9-b605-3113b9268794
# ╠═6e1ea2d0-6dd2-47ce-893d-1da75f14abc1
# ╠═6d2c0ec7-4b4a-409a-98e5-dc16c362558a
# ╠═8d9a068b-42aa-441c-86f8-7beab274f27c
# ╠═49546fa4-274d-4a62-8138-05046bfc7450
# ╠═d86e29cb-b67d-41b5-bb4a-9199a1e8b3db
# ╠═e45a2654-1a19-49c4-8a44-7b859dc0f5d1
# ╠═0923a7e7-4838-4d6a-babf-6551c1b0fa5d
# ╠═9d3afdde-7903-4539-ac3e-b56c38c4c834
# ╠═f539ff24-56cf-473d-a29f-66fafa579d50
# ╠═abbe0b94-f93c-4614-a534-020fc1526301
