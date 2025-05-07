using GLMakie
using Distributions
using LinearAlgebra

α = 0.1
δ = 0.1
γ = 0.01
β = 1e-3
η = 1.0
A = 1.0

# Compute C and S
C = (α * β * γ) / (β^2 + δ^2 * γ^2) - (α * β) / (β^2 + δ^2)
S = (α * δ * γ^2) / (β^2 + δ^2 * γ^2) - (α * δ) / (β^2 + δ^2)

# Compute sigma (standard deviation of noise)
σ² = (η^2 * α^2 * (1 - γ)^2) / (2 * δ * (1 + γ))
σ = sqrt(σ²)

# Time vector (extended to show more oscillations)
t_range = range(0, 2π/β, length=10000)

# Input function u(t) = A sin(βt)
u_t = A * sin.(β * t_range)

# Compute μ(t) = C cos(β t) - S sin(β t)
μ_t = C * cos.(β * t_range) - S * sin.(β * t_range)

# Compute P(t)
P_t = zeros(length(t_range))
for (i, t) in enumerate(t_range)
    if sign(μ_t[i]) == sign(cos(β * t))
        P_t[i] = cdf(Normal(), abs(μ_t[i]) / σ)
    else
        P_t[i] = cdf(Normal(), -abs(μ_t[i]) / σ)
    end
end

fig = with_theme(theme_latexfonts(), fontsize=24) do
    fig = Figure(size=(1200, 1000))

    ax_style = (;
        aspect=2,
        xtickalign=1,
        ytickalign=1,
    )

    # Create custom ticks at multiples of π/2
    custom_ticks = [k*π/2 for k in 0:8]
    tick_labels = ["0", "π/2", "π", "3π/2", "2π", "5π/2", "3π", "7π/2", "4π"]

    # Top plot: Input u(t)
    ax0 = Axis(fig[1, 1];
        xticklabelsvisible=false,
        ylabel=L"u(t)",
        ax_style...)
    lines!(ax0, β * t_range, u_t, color=:black, linewidth=4)
    xlims!(ax0, 0, β * maximum(t_range))
    ylims!(ax0, -1.1A, 1.1A)
    ax0.xticks = (custom_ticks, tick_labels)

    # Middle plot: μ(t) with uncertainty bands
    ax1 = Axis(fig[2, 1];
        xticklabelsvisible=false,
        ylabel=L"\mu(t) \pm k \, \sigma",
        ax_style...)
    
    hlines!(ax1, [0.0]; color="#ff6666", linewidth=3, linestyle = :dash)
    band!(ax1, β * t_range, μ_t .- 3σ, μ_t .+ 3σ, color=(:gray, 0.25))
    band!(ax1, β * t_range, μ_t .- σ, μ_t .+ σ, color=(:gray, 0.5))
    lines!(ax1, β * t_range, μ_t, color=:black, linewidth=4)
    
    xlims!(ax1, 0, β * maximum(t_range))
    ax1.xticks = (custom_ticks, tick_labels)

    # Bottom plot: P(t)
    ax2 = Axis(fig[3, 1];
        xlabel=L"\beta \, t",
        ylabel=L"\mathbb{P}(t)",
        ax_style...)
    
    lines!(ax2, β * t_range, P_t, color=:black, linewidth=4)
    hlines!(ax2, [0.5], color=:red, linestyle=:dash, linewidth=3)
    
    xlims!(ax2, 0, β * maximum(t_range))
    ylims!(ax2, -0.025, 1.025)
    ax2.xticks = (custom_ticks, tick_labels)

    # Adjust layout
    rowgap!(fig.layout, 1, 25)
    rowgap!(fig.layout, 2, 25)
    fig
end
