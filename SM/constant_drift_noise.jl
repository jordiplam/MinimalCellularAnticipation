using GLMakie
using Distributions

α = 0.1
δ = 0.1
γ = 0.01
β = 0.1
u0 = 0.0
x0 = 0.0
y0 = 0.0
η = 1.0

t_range = 0:0.1:25

function μ(t)
    term1 = (α * β) / (δ^2 * γ) * (1 - γ)
    term2 = (α * β / δ^2 - α * u0 / δ + x0) * exp(-δ * t)
    term3 = (α * β / (δ^2 * γ) - α * u0 / δ + y0) * exp(-δ * γ * t)
    return term1 + term2 - term3
end

function σ²(t)
    term1 = (η^2 * α^2) / (2 * δ) * (1 - exp(-2 * δ * t))
    term2 = (η^2 * γ * α^2) / (2 * δ) * (1 - exp(-2 * γ * δ * t))
    term3 = (2 * η^2 * γ * α^2) / (δ * (1 + γ)) * (1 - exp(-δ * (1 + γ) * t))
    return term1 + term2 - term3
end

# Calculate μ and σ for all time points
mu_vals = μ.(t_range)
sigma_vals = sqrt.(σ².(t_range))

# Create input function u(t)
u(t) = β * t + u0
u_vals = u.(t_range)

fig = with_theme(theme_latexfonts(), fontsize = 24) do
    fig = Figure(size = (1200, 400))

    ax_style = (;
        aspect = 2,
        xtickalign = 1,
        ytickalign = 1,
    )

    # Top plot: input u(t)
    ax1 = Axis(fig[1, 1];
        #xticklabelsvisible = false,
        ylabel = L"u(t)",
        ax_style...)
    lines!(ax1, t_range, u_vals, color=:black, linewidth=4)
    xlims!(ax1, extrema(t_range)...)

    # Bottom plot: μ(t) with uncertainty bands
    ax2 = Axis(fig[1, 2];
        xlabel = L"\text{Time } (t)",
        ylabel = L"\mu(t) \pm k \, \sigma(t)",
        ax_style...)
    
    hlines!(ax2, [0.0]; color="#ff6666", linewidth=3, linestyle = :dash)
    
    # Plot 3STD band
    band!(ax2, t_range, 
        mu_vals - 3*sigma_vals, 
        mu_vals + 3*sigma_vals, 
        color=(:gray, 0.25))
    
    # Plot 1STD band
    band!(ax2, t_range, 
        mu_vals - sigma_vals, 
        mu_vals + sigma_vals, 
        color=(:gray, 0.5))
    
    # Plot mean
    lines!(ax2, t_range, mu_vals, color=:black, linewidth=4)
    
    xlims!(ax2, extrema(t_range)...)
    ylims!(ax2, -0.5, 2.5)
    
    fig
end