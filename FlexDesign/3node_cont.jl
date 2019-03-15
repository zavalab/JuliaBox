using JuMP, Distributions, Gurobi, JSON

# Set the covariance matrix for the uncertain parameters
θ_nom = [0.; 60.; 10.]
β = 0.
covar = [80. β β; β 80. β; β β 120.]

# Set the dimensions
n_z = 3
n_θ = 3
n_d = 3

# Set the problem parameters
c = ones(n_d) / sqrt(n_d)
c_max = 0:0.25:10.5
U = 10000
num_samples = 1000
sig_dig = 6

# Obtain samples
srand(11)
d = MvNormal(θ_nom, covar)
θs = rand(d, num_samples)

# Setup storage arrays
cont_costs = zeros(length(c_max))
cont_SFs = zeros(length(c_max))
cont_times = zeros(length(c_max))
cont_objs = zeros(length(c_max))
feasibles = Vector{Bool}(length(c_max))
cont_data = Dict()

# Iterate and solve over the different costs
for itr = 1:length(c_max)
    # Initialize the main model and the variables
    m = Model(solver = GurobiSolver(OutputFlag = 0))
    @variable(m, 0 <= y[1:num_samples] <= 1)
    @variable(m, z[1:n_z, 1:num_samples])
    @variable(m, d[1:n_d] >= 0)

    # Set objective function
    @objective(m, Max, 1 / num_samples * sum(1 - y[k] for k = 1:num_samples))

    # Set the line capacity constraints
    @constraint(m, f1[k = 1:num_samples], -z[1, k] - 35 - d[1] <= y[k] * U)
    @constraint(m, f2[k = 1:num_samples], z[1, k] - 35 - d[1] <= y[k] * U)
    @constraint(m, f3[k = 1:num_samples], -z[2, k] - 50 - d[2] <= y[k] * U)
    @constraint(m, f4[k = 1:num_samples], z[1, k] - 50 - d[2] <= y[k] * U)

    # Set the generator capacity constraints
    @constraint(m, f5[k = 1:num_samples], -z[3, k] <= y[k] * U)
    @constraint(m, f6[k = 1:num_samples], z[3, k] - 100 - d[3] <= y[k] * U)

    # Set the node balance constraints
    @constraint(m, h1[k = 1:num_samples], z[1, k] - θs[1, k] == 0)
    @constraint(m, h2[k = 1:num_samples], -z[1, k] -z[2, k] + z[3, k] - θs[2, k] == 0)
    @constraint(m, h3[k = 1:num_samples], z[2, k] - θs[3, k] == 0)

    # Enforce the minimum SF
    @constraint(m, max_cost, sum(c[i] * d[i] for i = 1:n_d) <= c_max[itr])

    # Solve and and obtain results
    status = solve(m)
    opt_y = getvalue(y)
    opt_d = getvalue(d)
    opt_obj = getobjectivevalue(m)
    opt_time = getsolvetime(m)

    # Estimate the value of SF
    SF = 1 - sum(opt_y .>= 1e-8) / num_samples

    # Print the results
    print("------------------RESULTS------------------\n")
    print("Optimal Objective:     ", signif(opt_obj, sig_dig), "\n")
    print("Solution Time:         ", signif(opt_time, sig_dig), "s\n")
    print("Optimal Cost:          ", signif(sum(c[i] * opt_d[i] for i = 1:n_d), sig_dig), "\n")
    print("Maximum Cost:          ", signif(c_max[itr], sig_dig), "\n")
    print("Predicted SF:          ", 100 * signif(SF, sig_dig), "%\n")
    print("Optimal Design Values: ", opt_d, "\n\n")

    # Check solution feasibility for mapped y's
    yf = zeros(num_samples)
    yf[opt_y .>= 1e-8] = 1
    m2 = Model(solver = GurobiSolver(OutputFlag = 0))
    @variable(m2, z2[1:n_z, 1:num_samples])
    @variable(m2, d2[1:n_d] >= 0)

    # Set objective function
    @objective(m2, Max, 1 / num_samples * sum(1 - yf[k] for k = 1:num_samples))

    # Set the line capacity constraints
    @constraint(m2, f1[k = 1:num_samples], -z2[1, k] - 35 - d2[1] <= yf[k] * U)
    @constraint(m2, f2[k = 1:num_samples], z2[1, k] - 35 - d2[1] <= yf[k] * U)
    @constraint(m2, f3[k = 1:num_samples], -z2[2, k] - 50 - d2[2] <= yf[k] * U)
    @constraint(m2, f4[k = 1:num_samples], z2[1, k] - 50 - d2[2] <= yf[k] * U)

    # Set the generator capacity constraints
    @constraint(m2, f5[k = 1:num_samples], -z2[3, k] <= yf[k] * U)
    @constraint(m2, f6[k = 1:num_samples], z2[3, k] - 100 - d2[3] <= yf[k] * U)

    # Set the node balance constraints
    @constraint(m2, h1[k = 1:num_samples], z2[1, k] - θs[1, k] == 0)
    @constraint(m2, h2[k = 1:num_samples], -z2[1, k] -z2[2, k] + z2[3, k] - θs[2, k] == 0)
    @constraint(m2, h3[k = 1:num_samples], z2[2, k] - θs[3, k] == 0)

    # Enforce the minimum SF
    @constraint(m2, sum(c[i] * d2[i] for i = 1:n_d) <= c_max[itr])

    # Solve and and obtain results
    status = solve(m2)

    # Save results
    feasibles[itr] = status == :Optimal
    cont_objs[itr] = opt_obj
    cont_costs[itr] = sum(c[i] * opt_d[i] for i = 1:n_d)
    cont_SFs[itr] = SF
    cont_times[itr] = opt_time
    cont_data[string(itr)] = Dict("y" => opt_y, "d" => opt_d, "cost" => cont_costs[itr], "SF" => SF,
                                  "time" => opt_time, "obj" => opt_obj, "feasible" => feasibles[itr])
end

# Write results to data file
open("3d_cont_data2.json", "w") do f
    JSON.print(f,cont_data)
end
