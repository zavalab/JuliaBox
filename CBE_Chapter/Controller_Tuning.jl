using InfiniteOpt, Distributions, Ipopt, Random, Plots, Statistics

# Set constants
α = 0.001 # prioritize setpoint tracking
Random.seed!(42) # make the results reproducible

# Define the model
m = InfiniteModel(Ipopt.Optimizer)

# Add the infinite parameters
@infinite_parameter(m, t ∈ [0, 10], num_supports = 101)
@infinite_parameter(m, ysp ~ Normal(0, 2), num_supports = 100)

# Add the infinite variables (states and inputs)
@variable(m, y, Infinite(t, ysp))
@variable(m, -5 <= u <= 5, Infinite(t, ysp))
@variable(m, e, Infinite(t, ysp))
@variable(m, I, Infinite(t, ysp))

# Add the finite variables (the PID tuning parameters)
@variable(m, 0 <= Kc <= 10)
@variable(m, 0 <= τI <= 10)
@variable(m, 0 <= τD <= 10)

# Add the objective
@objective(m, Min, 𝔼(∫(e^2 + α * u^2, t), ysp))

# Add the constraints
@constraint(m, ∂(y, t) == -y + u)
@constraint(m, u == Kc * e + τI * I - τD * ∂(y, t))
@constraint(m, e == ysp - y)
@constraint(m, ∂(I, t) == e)
@constraint(m, y(0, ysp) == 0)
@constraint(m, I(0, ysp) == 0)

# Optimize and get the results
optimize!(m)
println("Kc = ", value(Kc))
println("τI = ", value(τI))
println("τD = ", value(τD))

# Plot the tracking error statistics
ts = value(t)
e_opt = value(e, ndarray = true)
e_mean = mean(e_opt, dims = 2)
e_var = var(e_opt, dims = 2)
plot(ts, e_mean, ribbon = e_var, fillalpha = 0.2, label = "Mean Error")
xlabel!("Time (t)")
ylabel!("Tracking Error (e(t))")
