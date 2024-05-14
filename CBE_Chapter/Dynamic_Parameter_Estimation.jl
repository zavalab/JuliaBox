using InfiniteOpt, Ipopt, Plots

# Define the experimental time series data
t_exp = collect(0:2.4:24)
y_exp = [0.0106   0.0120  0.0170  0.0202  0.025  0.031  0.039  0.058  0.077  0.083  0.111
         0.0091   0.0147  0.0245  0.0358  0.045  0.081  0.092  0.133  0.153  0.160  0.168]

# Create the InfiniteModel objective
m = InfiniteModel(Ipopt.Optimizer)

# Define the time domain for simulation
# adds 4 (i.e., 6 - 2) collocation pts between supports to approximate dydt
dmethod = OrthogonalCollocation(6)
@infinite_parameter(m, t ∈ [0, 24], supports = t_exp, derivative_method = dmethod)

# Define variables
@variable(m,  0 <= r[i in 1:2] <= 10, start = 1) # growth rates
@variable(m, -3 <= a[i in 1:2, k in 1:2] <= 3, start = 1) #  interaction parameters
@variable(m,  0 <= y[i in 1:2] <= 1, Infinite(t), start = 0.5) # abundance

# Define the model
@constraint(m, [i in 1:2], ∂(y[i], t) == (r[i] + sum(a[i, j] * y[j] for j in 1:2)) * y[i])

# Define the objective
@objective(m, Min,
    sum((y[i](t_exp[k]) - y_exp[i, k])^2 for i in 1:2, k in 1:length(t_exp))
)

# Solve and get the results
optimize!(m)
r_opt = value.(r)
a_opt = value.(a)
t_opt = value(t, label = All) # extract all the extra time points that were used
y_opt = value.(y, label = All)

# Display the results
println("r_1 = ", r_opt[1], ", a_11 = ", a_opt[1, 1], ", a_12 = ", a_opt[1, 2])
println("r_2 = ", r_opt[2], ", a_21 = ", a_opt[2, 1], ", a_22 = ", a_opt[2, 2])
plot(t_opt, y_opt, label = ["Species 1" "Species 2"],
    color = [:gray :black], linestyle = [:dash :solid], linewidth = [2 2],  legend = :best,
)
scatter!(t_exp, y_exp', label = ["Species 1 (exp)" "Species 2 (exp)"],
    color = [:gray :black],
    markershape = [:utriangle :circle]
)
xlabel!("Time (hr)")
ylabel!("Abundance")
