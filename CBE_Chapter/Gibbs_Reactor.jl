# Load in necessary packages
using JuMP, Ipopt

# Define JuMP Model Object
m = Model(Ipopt.Optimizer)

# Define set
K = ["CO", "H2", "CH3OH"]

# Define problem data
μ_in = [100, 600, 0] # Inlet molar flows (kmol/hr)
μ_in = Dict(zip(K, μ_in)) # Create dictionary with keys K and values μ_in
P = 150 # Pressure (bar)
T = 300 + 273.15 # Temperature (K)
Keq = 10^(-12.275+4938/T) # Equilibrium constant
γ = [-1, -2, 1] # Stochiometric Coefficients
γ = Dict(zip(K, γ)) # Create dictionary with keys K and values γ

# Define variables
@variable(m, μ_out[K] >= 0)
@variable(m, μ_tot >= 0)
@variable(m, a[K] >= 0)
@variable(m, ξ >= 0)

# Define constraints
@constraint(m, [k in K], μ_out[k] == μ_in[k] + γ[k] * ξ) # Conservation
@constraint(m, μ_tot == sum(μ_out[k] for k in K)) # Total flow
@constraint(m, [k in K], a[k] == (μ_out[k] / μ_tot)^ γ[k]) # Activities
@constraint(m, a["CO"] * a["H2"] * a["CH3OH"] == (P^2) * Keq) # Equilibrium condition

# Solve and display solution
optimize!(m)
println("ξ = ", value(ξ))
println("μ_out = ", value.( μ_out))
