using GLMakie

# Set parameters
α = 0.1
δ = 0.1
γ = 0.01

# Define the system response for piecewise linear input
function system_response_rate(β1, β2, ts, u0, t_range)
    u = zeros(length(t_range))
    x = zeros(length(t_range))
    y = zeros(length(t_range))
    
    for (i, t) in enumerate(t_range)
        # Input function
        if t < ts
            u[i] = β1 * t + u0
        else
            τ = t - ts
            u[i] = β2 * τ + β1 * ts + u0
        end
        
        # System response (approximate solution for t ≥ ts)
        if t < ts
            # Before switch - assume transient dynamics have settled
            x[i] = (α * (β1 * t + u0)) / δ - (α * β1) / δ^2
            y[i] = (α * (β1 * t + u0)) / δ - (α * β1) / (δ^2 * γ)
        else
            τ = t - ts
            # After switch - use the given approximate solution
            x[i] = (α * (β2 * (τ + ts) + (β1 - β2) * ts + u0)) / δ - (α * β2) / δ^2 +
                   (α * (β2 - β1)) / δ^2 * exp(-δ * τ)
            y[i] = (α * (β2 * (τ + ts) + (β1 - β2) * ts + u0)) / δ - (α * β2) / (δ^2 * γ) +
                   (α * (β2 - β1)) / (δ^2 * γ) * exp(-δ * γ * τ)
        end
    end
    
    return u, x, y
end

t_range = -500:1:1500
ts = 0

# Case 1: β1 < β2
β1_a = -0.01
β2_a = 0.01
u0_a = 1
u_a, x_a, y_a = system_response_rate(β1_a, β2_a, ts, u0_a, t_range)

# Case 2: β1 > β2
β1_b = 0.01
β2_b = -0.01
u0_b = 1
u_b, x_b, y_b = system_response_rate(β1_b, β2_b, ts, u0_b, t_range)

fig = with_theme(theme_latexfonts(), fontsize = 24) do
    fig = Figure(size = (1200, 650))

    ax_style = (;
        aspect = 2,
        xtickalign = 1,
        ytickalign = 1,
    )

    # Top left: u for β1 < β2
    ax1 = Axis(fig[1, 1];
        xticklabelsvisible = false,
        ylabel = L"u(\tau)",
        ax_style...)
    lines!(ax1, t_range, u_a, color=:black, linewidth=4)
    xlims!(ax1, extrema(t_range)...)

    # Bottom left: x-y for β1 < β2
    ax2 = Axis(fig[2, 1];
        xlabel = L"\tau",
        ylabel = L"x(\tau) - y(\tau)",
        ax_style...)
    hlines!(ax2, [0], color="#ff6666", linewidth=2)
    lines!(ax2, t_range, x_a - y_a, color=:gray, linestyle=(:dash, :dense), linewidth=4)
    xlims!(ax2, extrema(t_range)...)
    ylims!(ax2, -10, 10)

    # Top right: u for β1 > β2
    ax3 = Axis(fig[1, 2];
        xticklabelsvisible = false,
        yticklabelsvisible = false,
        ax_style...)
    lines!(ax3, t_range, u_b, color=:black, linewidth=4)
    xlims!(ax3, extrema(t_range)...)

    # Bottom right: x-y for β1 > β2
    ax4 = Axis(fig[2, 2];
        xlabel = L"\text{Time } (t)",
        yticklabelsvisible = false,
        ax_style...)
    hlines!(ax4, [0], color="#ff6666", linewidth=2)
    lines!(ax4, t_range, x_b - y_b, color=:gray, linestyle=(:dash, :dense), linewidth=4)
    xlims!(ax4, extrema(t_range)...)
    ylims!(ax4, -10, 10)

    # Adjust layout
    colgap!(fig.layout, 1, 35)
    rowgap!(fig.layout, 1, 5)
    fig
end