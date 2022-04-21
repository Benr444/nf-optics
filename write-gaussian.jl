using Images, Plots

edge_length = 768
function gaussian(x, y, x_0, y_0, σ_x, σ_y)
	diff_x = ((x - x_0) / σ_x)^2
	diff_y = ((y - y_0) / σ_y)^2
	return exp(-0.5 * (diff_x + diff_y))
end
function gaussian(x, μ, σ)
	diff = (x - μ / σ)^2
	scale = 1 / (σ * sqrt(2π))
	return scale * exp(-0.5 * diff)
end

X = zeros(edge_length, edge_length)
mean_x = mean_y = edge_length / 2.0
σ_x = σ_y = edge_length / 4.0 # empirical rule: 95%

for i in 1:edge_length
	for j in 1:edge_length
		X[i, j] = gaussian(i, j, mean_x, mean_y, σ_x, σ_y)
	end
end
	
#display(heatmap(X))

save("in/gaussian-$(edge_length)px.png", X)