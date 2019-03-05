using JuMP, Distributions, CPLEX, JSON

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
num_samples = 2000
sig_dig = 6

# Obtain samples
srand(11)
d = MvNormal(θ_nom, covar)
θs = rand(d, num_samples)

# Setup storage arrays
bin_costs = zeros(length(c_max))
bin_SFs = zeros(length(c_max))
bin_times = zeros(length(c_max))
bin_data = Dict()

# Setup the Model
m = Model(solver = CplexSolver(CPX_PARAM_TILIM = 3600))
@variable(m, y[1:num_samples], Bin)
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
@constraint(m, max_cost, sum(c[i] * d[i] for i = 1:n_d) <= 0)

# Iterate and solve over the different costs
for itr = 1:length(c_max)
    # Set the cost
    JuMP.setRHS(max_cost, c_max[itr])

    # Solve and and obtain results
    status = solve(m)
    opt_y = getvalue(y)
    opt_d = getvalue(d)
    opt_obj = getobjectivevalue(m)
    opt_time = getsolvetime(m)

    # Estimate the value of SF
    SF = opt_obj

    # Print the results
    print("------------------RESULTS------------------\n")
    print("Optimal Objective:     ", signif(opt_obj, sig_dig), "\n")
    print("Solution Time:         ", signif(opt_time, sig_dig), "s\n")
    print("Optimal Cost:          ", signif(sum(c[i] * opt_d[i] for i = 1:n_d), sig_dig), "\n")
    print("Maximum Cost:          ", signif(c_max[itr], sig_dig), "\n")
    print("Predicted SF:          ", 100 * signif(SF, sig_dig), "%\n")
    print("Optimal Design Values: ", opt_d, "\n\n")

    # Save results
    bin_costs[itr] = sum(c[i] * opt_d[i] for i = 1:n_d)
    bin_SFs[itr] = SF
    bin_times[itr] = opt_time
    bin_data[itr] = Dict("y" => opt_y, "d" => opt_d, "cost" => bin_costs[itr], "SF" => SF, "time" => opt_time)
end

# Write results to data file
open("ieee14_bin_data.json", "w") do f
    JSON.print(f,bin_data)
end
