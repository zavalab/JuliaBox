using JuMP, Distributions, Gurobi, JSON

# Set the covariance matrix for the uncertain parameters
θ_nom = [87.3; 50.0; 25.0; 28.8; 50.0; 25.0; 0; 0; 0; 0; 0]
β = 240.
covar = eye(length(θ_nom)) * 1200.
covar[covar .== 0] = β

# Specify the network details
line_cap = 100
gen_cap = [332; 140; 100; 100; 100]

# Set the dimensions
n_z = 25
n_θ = 11
n_d = 25

# Set the problem parameters
c = ones(n_d) / sqrt(n_d)
c_max = 0:2.5:65
U = 10000
num_samples = 2000 # Change to desired amount
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

    # Initialize the model and the variables
    m = Model(solver = GurobiSolver(OutputFlag = 0))
    @variable(m, 0 <= y[1:num_samples] <= 1)
    @variable(m, z[1:n_z, 1:num_samples])
    @variable(m, d[1:n_d] >= 0)

    # Set objective function
    @objective(m, Max, 1 / num_samples * sum(1 - y[k] for k = 1:num_samples))

    # Set the line capacity constraints
    @constraint(m, [line = 1:20, k = 1:num_samples], -z[line, k] - line_cap - d[line] <= y[k] * U)
    @constraint(m, [line = 1:20, k = 1:num_samples], z[line, k] - line_cap - d[line] <= y[k] * U)

    # Set the generator capacity constraints
    @constraint(m, [gen = 1:5, k = 1:num_samples], -z[gen + 20, k] <= y[k] * U)
    @constraint(m, [gen = 1:5, k = 1:num_samples], z[gen + 20, k] - gen_cap[gen] - d[gen + 20] <= y[k] * U)

    # Set the node balance constraints
    @constraint(m, h1[k = 1:num_samples], z[21, k] - z[1, k] - z[6, k] == 0)
    @constraint(m, h2[k = 1:num_samples], sum(z[i, k] for i = [22; 1]) - sum(z[i, k] for i = [2; 4; 5]) - θs[1, k] == 0)
    @constraint(m, h3[k = 1:num_samples], sum(z[i, k] for i = [23; 2]) - z[3, k] - θs[2, k] == 0)
    @constraint(m, h4[k = 1:num_samples], sum(z[i, k] for i = [3; 4; 8]) - sum(z[i, k] for i = [7; 11]) - θs[3, k] == 0)
    @constraint(m, h5[k = 1:num_samples], sum(z[i, k] for i = [5; 6; 7; 12]) - θs[4, k] == 0)
    @constraint(m, h6[k = 1:num_samples], sum(z[i, k] for i = [24; 16; 18]) - sum(z[i, k] for i = [12; 19]) - θs[5, k] == 0)
    @constraint(m, h7[k = 1:num_samples], z[9, k] - sum(z[i, k] for i = [8; 10]) == 0)
    @constraint(m, h8[k = 1:num_samples], z[25, k] - z[9, k] == 0)
    @constraint(m, h9[k = 1:num_samples], sum(z[i, k] for i = [10; 11]) - sum(z[i, k] for i = [13; 14]) - θs[6, k] == 0)
    @constraint(m, h10[k = 1:num_samples], sum(z[i, k] for i = [13; 20]) - θs[7, k] == 0)
    @constraint(m, h11[k = 1:num_samples], z[19, k] - z[20, k] - θs[8, k] == 0)
    @constraint(m, h12[k = 1:num_samples], z[17, k] - z[18, k] - θs[9, k] == 0)
    @constraint(m, h13[k = 1:num_samples], z[15, k] - sum(z[i, k] for i = [16; 17]) - θs[10, k] == 0)
    @constraint(m, h14[k = 1:num_samples], z[14, k] - z[15, k] - θs[11, k] == 0)

    # Enforce the max cost
    @constraint(m, max_cost, sum(c[i] * d[i] for i = 1:n_d) <= c_max[itr])

    # Solve and and obtain results
    status = solve(m)
    opt_y = getvalue(y)
    opt_d = getvalue(d)
    opt_obj = getobjectivevalue(m)
    opt_time = getsolvetime(m)

    # Estimate the value of SF
    SF = 1 - sum(opt_y .>= 1e-5) / num_samples

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
    yf[opt_y .>= 1e-6] = 1
    m2 = Model(solver = GurobiSolver(OutputFlag = 0))
    @variable(m2, z2[1:n_z, 1:num_samples])
    @variable(m2, d2[1:n_d] >= 0)

    # Set objective function
    @objective(m2, Max, 1 / num_samples * sum(1 - yf[k] for k = 1:num_samples))

    # Set the line capacity constraints
    @constraint(m2, [line = 1:20, k = 1:num_samples], -z2[line, k] - line_cap - d2[line] <= yf[k] * U)
    @constraint(m2, [line = 1:20, k = 1:num_samples], z2[line, k] - line_cap - d2[line] <= yf[k] * U)

    # Set the generator capacity constraints
    @constraint(m2, [gen = 1:5, k = 1:num_samples], -z2[gen + 20, k] <= yf[k] * U)
    @constraint(m2, [gen = 1:5, k = 1:num_samples], z2[gen + 20, k] - gen_cap[gen] - d2[gen + 20] <= yf[k] * U)

    # Set the node balance constraints
    @constraint(m2, h1[k = 1:num_samples], z2[21, k] - z2[1, k] - z2[6, k] == 0)
    @constraint(m2, h2[k = 1:num_samples], sum(z2[i, k] for i = [22; 1]) - sum(z2[i, k] for i = [2; 4; 5]) - θs[1, k] == 0)
    @constraint(m2, h3[k = 1:num_samples], sum(z2[i, k] for i = [23; 2]) - z2[3, k] - θs[2, k] == 0)
    @constraint(m2, h4[k = 1:num_samples], sum(z2[i, k] for i = [3; 4; 8]) - sum(z2[i, k] for i = [7; 11]) - θs[3, k] == 0)
    @constraint(m2, h5[k = 1:num_samples], sum(z2[i, k] for i = [5; 6; 7; 12]) - θs[4, k] == 0)
    @constraint(m2, h6[k = 1:num_samples], sum(z2[i, k] for i = [24; 16; 18]) - sum(z2[i, k] for i = [12; 19]) - θs[5, k] == 0)
    @constraint(m2, h7[k = 1:num_samples], z2[9, k] - sum(z2[i, k] for i = [8; 10]) == 0)
    @constraint(m2, h8[k = 1:num_samples], z2[25, k] - z2[9, k] == 0)
    @constraint(m2, h9[k = 1:num_samples], sum(z2[i, k] for i = [10; 11]) - sum(z2[i, k] for i = [13; 14]) - θs[6, k] == 0)
    @constraint(m2, h10[k = 1:num_samples], sum(z2[i, k] for i = [13; 20]) - θs[7, k] == 0)
    @constraint(m2, h11[k = 1:num_samples], z2[19, k] - z2[20, k] - θs[8, k] == 0)
    @constraint(m2, h12[k = 1:num_samples], z2[17, k] - z2[18, k] - θs[9, k] == 0)
    @constraint(m2, h13[k = 1:num_samples], z2[15, k] - sum(z2[i, k] for i = [16; 17]) - θs[10, k] == 0)
    @constraint(m2, h14[k = 1:num_samples], z2[14, k] - z2[15, k] - θs[11, k] == 0)

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
open("ieee14_cont_data2.json", "w") do f
    JSON.print(f,cont_data)
end
