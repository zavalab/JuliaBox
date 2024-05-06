using InfiniteOpt, Distributions, LinearAlgebra, KNITRO

# Define constants
α, β,  M, tf = 0.95, 1.3, 100, 10
μ = zeros(10)
Σ = I * 3

# Initialize the model
m = InfiniteModel(KNITRO.Optimizer)

# Add the infinite parameters corresponding to the infinite domains
@infinite_parameter(m, t ∈ [0, tf], num_supports = 10)
@infinite_parameter(m, ξ[1:10] ~ MvNormal(μ, Σ), num_supports = 10)

# Add the variables and their domain constraints
@variable(m, ya ≥ 0, Infinite(t))
@variable(m, yb ≥ 0, Infinite(t, ξ))
@variable(m, yc, Bin, Infinite(ξ))
@variable(m, z[1:2], Int)

# Define the objective
@objective(m, Min, ∫(ya^2 + 2 * 𝔼(yb, ξ), t))

# Add the constraints
@constraint(m, ∂(yb, t) == 2yb^2 + ya - z[1])
@constraint(m, yb ≤ yc * M)
@constraint(m, 𝔼(yc, ξ) ≥ α)
@constraint(m, ya(0) + z[2] == β)

# Solve and retrieve the results
optimize!(m)
println("Objective value: ", objective_value(m))
println("ya(t): ", value(ya))
println("z: ", value.(z))
