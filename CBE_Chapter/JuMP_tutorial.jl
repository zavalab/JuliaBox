# Load in necessary packages
using JuMP, Ipopt

# Define JuMP Model Object
m = Model(Ipopt.Optimizer)
set_optimizer_attribute(m, "max_iter", 100)
set_optimizer_attribute(m, "tol", 1e-8)

# Define parameter data
p = 2

# Define variables
@variable(m, x1 >= 0)
@variable(m, x2 >= 0, start = 1)

# Define constraints
@constraint(m, x1 + x2 == p)
@constraint(m, x1 / x2 == 1)

# Define objective
@objective(m, Min, x1 + x2^2)

# Solve model
optimize!(m)

# Display results
println("Status = ", termination_status(m))
println("x1 = ", value(x1))
println("x2 = ", value(x2))
println("ϕ = ", objective_value(m))



### Set Notation Example

# Load in necessary packages
using JuMP, Ipopt

# Define JuMP Model Object
m = Model(Ipopt.Optimizer)

# Define set
K = ["myvar1", "myvar2"]

# Define variables
@variable(m, x[K] >= 0, start = 1)

# Define constraints
@constraint(m, sum(x[k] for k in K) == 2)
@constraint(m, x["myvar1"] / x["myvar2"] == 1)

# Define objective
@objective(m, Min, x["myvar1"] + x["myvar2"]^2)

# Solve model
optimize!(m)

# Display results
println("Status = ", termination_status(m))
println("x = ", value.(x))
println("ϕ = ", objective_value(m))

### Sensitivity Analysis
using JuMP, Ipopt, PyPlot

# Define the JuMP model with a parameter object
optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 100, "tol" => 1e-8)
m = Model(optimizer)

@variable(m, x1 >= 0)
@variable(m, x2 >= 0, start = 1)
@variable(m, p in Parameter(1))

@constraint(m, x1 + x2 == p)
@constraint(m, x1 / x2 == 1)

@objective(m, Min, x1 + x2^2)

# Create vector spanning parameter range and placeholder for optimal objective
ps = 1:1:10
ϕs = zeros(length(ps))

# Run loop to solve problem
for j = 1:length(ps)
    set_parameter_value(p, ps[j])
    optimize!(m)
    ϕs[j] = objective_value(m)
end

plot(ps, ϕs, color = "blue"); grid("on"); xlabel("p"); ylabel(L"$\varphi$")
savefig((@__DIR__)*"/plot.pdf")
