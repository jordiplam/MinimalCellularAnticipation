using DataInterpolations
using DiffEqNoiseProcess
using GLMakie
using ModelingToolkit
using OrdinaryDiffEq
using RandomNumbers

# Import functions from ModelingToolkit with clearer names.
using ModelingToolkit: t_nounits as t, D_nounits as D

"""
	generate_u(tspan, W_dt, λ; seed=0, saveat=nothing)

Generate a smoothed time series based on a Wiener process.

This function creates a Wiener process sampled at intervals W_dt, then applies
exponential smoothing with parameter λ to create a continuous, non-negative
function u(t).

# Arguments
- `tspan`: Tuple (t_start, t_end) for the simulation time range.
- `W_dt`: Time step for Wiener process calculation.
- `λ`: Smoothing parameter (0 ≤ λ ≤ 1). Higher values give more smoothing.

# Kwargs
- `seed`: Random seed for reproducibility (default: 0).
- `saveat`: Optional time points to save (default: nothing).

# Returns
A LinearInterpolation object representing the smoothed function u(t).
"""
function generate_u(tspan, W_dt, λ; seed = 0, saveat = nothing)
	# Initialize Wiener process with specified random seed.
	W = WienerProcess(first(tspan), 0.0, nothing;
		rng = Xorshifts.Xoroshiro128Plus(seed),
	)

	# Create a realization of the Wiener process.
	# May overcompute but it is better to lack the last datapoint.
	while W.t[end] <= last(tspan)
		calculate_step!(W, W_dt, nothing, nothing)
		accept_step!(W, W_dt, nothing, nothing)
	end

	# Get time and shift the input to avoid negative concentrations.
	t = W.t
	u = W.u .- minimum(W.u)

	# Apply exponential smoothing.
	for i = 2:lastindex(u)
		u[i] = λ*u[i] + (1 - λ)*u[i - 1]
	end

	# Return as continuous interpolation.
	LinearInterpolation(u, t)
end

"""
	solve_system(sys, u0, tspan, ps, W_dt, λ; seed=0, saveat=[])

Solve the ODE system with the generated u(t) as an input.

# Arguments
- `sys`: ModelingToolkit ODESystem.
- `u0`: Initial conditions.
- `tspan`: Simulation time range.
- `ps`: System parameters.
- `W_dt`: Wiener process time step.
- `λ`: Smoothing parameter.

# Kwargs
- `seed`: Random seed (default: 0).
- `saveat`: Time points to save solution (default: []).

# Returns
ODESolution object.
"""
function solve_system(sys, u0, tspan, ps, W_dt, λ; seed = 0, saveat = [])
	# Create ODE problem with seeded random input u(t) as a parameter.
	prob = ODEProblem(sys, u0, tspan, [
			ps...,
			:u_t => generate_u(tspan, W_dt, λ; seed)
		],
	)

	# Solve problem.
	solve(prob, AutoTsit5(Rosenbrock23()); saveat)
end

"""
	trend(u)

Calculate the trend (direction of change) of a time series.

# Arguments
- `u`: Input time series

# Returns
Array of +1 (increasing) or -1 (decreasing) values
"""
function trend(u)
	sign.(diff(u))
end

"""
	predicted_trend(P)

Convert probability values to predicted trends.

# Arguments
- `P`: Probability values (0 ≤ P ≤ 1)

# Returns
Array of +1 (P > 0.5) or -1 (P < 0.5) values
"""
function predicted_trend(P)
	sign.(P .- 0.5)
end

"""
	main()

Main function to execute the full analysis pipeline.

- Defines the ODE system
- Solves the system with generated noise input
- Analyzes trends
- Creates visualization

# Returns
Figure object showing the results
"""
function main()
	# Define symbolic variables and parameters.
	@variables u(t) x(t) y(t) P(t)
	@parameters α γ δ (u_t::LinearInterpolation)(..)
	
	# Define system equations.
	eqs = [
		P	~ x/(x + y),
		u	~ u_t(t),
		D(x) ~ α*u - δ*x,
		D(y) ~ γ*(α*u - δ*y),
	]
	
	# Build the ODE system
	@mtkbuild sys = ODESystem(eqs, t)

	# Simulation parameters
	tspan = (-1e-3, 1e4)
	u0 = [u => 0 for u in unknowns(sys)]
	ps = [
		:α	=> 1e-1,
		:γ	=> 1e-2,
		:δ	=> 1e-1,
	]

	# Noise generation parameters.
	W_dt = 1e1
	λ = 2e-1

	# Solve the system.
	sol = solve_system(sys, u0, tspan, ps, W_dt, λ;
		seed = 0,
		saveat = 0:1e-2:last(tspan),
	)

	# Extract concentrations and trends.
	u_vals = sol(sol.t; idxs = sys.u).u
	xtrend = trend(u_vals)
	P_vals = sol(sol.t; idxs = sys.P).u
	ytrend = predicted_trend(P_vals)

	# Intitialize figure.
	fig = Figure()
	ax = Axis(fig[1, 1];
		aspect = 2,
		xlabel = "Time (t)",
		ylabel = "Concentration",
	)
	
	# Find regions where predicted trend changes.
	trend_change_idx = vcat(
		[1],
		findall(x -> x != 0, diff(ytrend)),
		length(ytrend) + 1,
	)
	trend_change_span = zip(
		trend_change_idx[1:end - 1],
		trend_change_idx[2:end] .- 1,
	)
	
	# Create colored regions showing predicted trends.
	trend_start = [sol.t[i] for (i, j) in trend_change_span]
	trend_end   = [sol.t[j] for (i, j) in trend_change_span]
	trend_value = [P_vals[i + 1] for (i, j) in trend_change_span]
	vspan!(trend_start, trend_end;
		color = ifelse.(trend_value .> 0.5, :blue, :red),
		alpha = 0.1,
	)
	
	# Plot input and predictions.
	lines!(ax, sol.t, u_vals/maximum(u_vals);
		label = L"\text{Normalized } u(t)",
		color = :black,
	)
	lines!(ax, sol.t, P_vals;
		label = L"P^*(t)",
		color = :gray,
	)
	xlims!(ax, tspan...)

	axislegend(ax)
	
	return fig
end

fig = main()