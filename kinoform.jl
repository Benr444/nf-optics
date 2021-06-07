### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ e2d9a9eb-b53d-45cd-95f9-2b2fe817480a
using Plots

# ╔═╡ 51b409ce-adc6-4dec-b8bf-419662784847
main_dir = "D:\\records\\usc\\magnetics-optics\\code\\r-optics\\"

# ╔═╡ 158dd0ce-c7ac-11eb-1e2c-9dd39891e973
include(main_dir * "ROptics.jl")

# ╔═╡ 88327064-4c98-480c-a825-0a927dad88ca
gr()

# ╔═╡ 3b8ea9a0-0597-4567-84c5-0030b69c00fb
default(
	xtickfont=30,
	ytickfont=30,
	size=(1500, 1500),
	fmt = :png
)

# ╔═╡ 82716bd6-f66d-466a-a575-cccf66c7d73d
my_λ = 1.0 # units:

# ╔═╡ f334efc6-c537-4359-b89d-9a6aa7e1df83
my_scale = 1.0 # units:

# ╔═╡ 914ef9ad-6e12-4401-84fd-e9f1639ed550
g, gˈ, G, Gˈ = Dict(), Dict(), Dict(), Dict() # prime marker is \verts

# ╔═╡ 89cb6f53-1a02-4e9e-9b3c-017873b8755c
G[0] = ROptics.OpticField(main_dir * "in/three.png")

# ╔═╡ 31b1ecf5-e7ad-4125-8f0e-d1a89e44caf5
g[0] = ROptics.F⁻¹(G[0])

# ╔═╡ d9b1e557-70b9-44ec-afb0-23e4736fe1ae
ROptics.png(g[0], main_dir * "out/three-g0.png")

# ╔═╡ c421bf2f-1912-4fc8-a936-bf38f6a33a85
total_cycles = 3

# ╔═╡ 0911ddb1-e8bf-4008-ad80-434c07374e1b
ROptics.OpticField

# ╔═╡ ecca5711-377d-4e62-b966-568fcb8c70d5
ROptics

# ╔═╡ Cell order:
# ╠═51b409ce-adc6-4dec-b8bf-419662784847
# ╠═158dd0ce-c7ac-11eb-1e2c-9dd39891e973
# ╠═e2d9a9eb-b53d-45cd-95f9-2b2fe817480a
# ╟─88327064-4c98-480c-a825-0a927dad88ca
# ╠═3b8ea9a0-0597-4567-84c5-0030b69c00fb
# ╠═82716bd6-f66d-466a-a575-cccf66c7d73d
# ╠═f334efc6-c537-4359-b89d-9a6aa7e1df83
# ╠═914ef9ad-6e12-4401-84fd-e9f1639ed550
# ╠═89cb6f53-1a02-4e9e-9b3c-017873b8755c
# ╠═31b1ecf5-e7ad-4125-8f0e-d1a89e44caf5
# ╠═d9b1e557-70b9-44ec-afb0-23e4736fe1ae
# ╠═c421bf2f-1912-4fc8-a936-bf38f6a33a85
# ╠═0911ddb1-e8bf-4008-ad80-434c07374e1b
# ╠═ecca5711-377d-4e62-b966-568fcb8c70d5
