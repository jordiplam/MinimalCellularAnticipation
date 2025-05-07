using GLMakie
using Roots

import Makie

function f(τ, α, δ, γ, β1, β2)
	factor = (α / (δ^2*γ))
	return factor*(β2*(1 - γ) + (β2 - β1)*(γ*exp(-δ * τ) - exp(-δ*γ*τ)))
end

function find_root(α, δ, γ, β1, β2)
	try
		result = find_zero(τ -> f(τ, α, δ, γ, β1, β2), 0)
		return result
	catch
		return NaN
	end
end

function create_heatmap_rootfinding(α, δ, γ, β1_range, β2_range)
	heatmap_data = zeros(length(β1_range), length(β2_range))
	for (i, β1) in enumerate(β1_range)
		for (j, β2) in enumerate(β2_range)
			heatmap_data[i, j] = find_root(α, δ, γ, β1, β2)
		end
	end
	heatmap_data
end

# Define parameters and ranges
α = 1e-1
δ = 1e-1
γ = 1e-2
β1_range = 10 .^ (range(log10(1e-4), log10(1), 500))
β1_range = unique(sort([-β1_range..., -1:0.01:1..., β1_range...]))
β2_range = 10 .^ (range(log10(1e-4), log10(1), 500))
β2_range = unique(sort([-β2_range..., -1:0.01:1..., β2_range...]))

# Create and display the heatmap
heatmap = create_heatmap_rootfinding(α, δ, γ, β1_range, β2_range)

fig = with_theme(theme_latexfonts(), fontsize = 30) do
	filtered_heatmap = copy(heatmap)
	filtered_heatmap[heatmap .>= 1e4] .= 1e4
	fig = Figure(size = (1000, 700))
	fticks = (
		-1:0.5:1,
		["-1", "-0.5", "0", "0.5", "1"],
	)
	ax = Axis3(fig[1, 1];
		xlabel = L"\beta_2",
		ylabel = L"\beta_1",
		zlabel = L"\tau_\text{r}",
		xlabeloffset = 75,
		ylabeloffset = 75,
		zlabeloffset = 100,
		zlabelrotation = 0,
		xticks = fticks,
		yticks = fticks,
		zticks = (2500:2500:10000, string.(2500:2500:10000)),
	)
	xlims!(ax, -1, 1)
	ylims!(ax, -1, 1)
	zlims!(ax, 0, 1e4)
	hm = surface!(ax, β2_range, β1_range, filtered_heatmap';
		colorscale = Makie.pseudolog10,
		rasterize = 10,
		shading = NoShading,
	)
	fig
end
