using FlexJuMP, JuMP, Gurobi, Pavito, Ipopt

include("data_reader.jl") # provides data reader function

# Set the dimensions
n_z = 141
n_θ = 84
n_d = 141
n_f = 282
n_h = 141

# Setup the parameters
covar = eye(n_θ) * 100

# Specify the network details
line_cap = 100
data = ParseFile("141bus.txt", line_cap)

# Extract information
fConsts = data["fConsts"]
fControls = data["fControls"]
fRandoms = data["fRandoms"]
hConsts = data["hConsts"]
hControls = data["hControls"]
hRandoms = data["hRandoms"]

# Setup the model
m = FlexibilityModel(solver = PavitoSolver(mip_solver = GurobiSolver(OutputFlag = 0),
                     cont_solver = IpoptSolver(print_level = 0), log_level = 0, mip_solver_drives = false))

# Define variables
@randomvariable(m, θ[i = 1:n_θ, mean = 0) # Temperarily set the mean to 0
@recoursevariable(m, z[1:n_f])

# Set the capacity constraints
@constraint(m, caps[j = 1:n_f], fConsts[j] + sum(fControls[j, i] * z[i] for i = 1:n_z) <= 0)

# Set the balances
@constraint(m, bal[j = 1:n_h], sum(hControls[j, i] * z[i] for i = 1:n_z) + sum(hRandoms[j, i] * θ[i] for i = 1:n_θ) == 0)

# Define the uncertainty set
setuncertaintyset(m, :Ellipsoid, covar)

# Compute the center
new_mean = findcenteredmean(m, center = :feasible, update_mean = true)

# Rank the constraints of design A and print
rank_data = rankinequalities(m, max_ranks = 300)
println("****Ranking Results****")
for i = 1:length(rank_data)
    println("Rank ", i)
    println("  Active Constraints: ", rank_data[i]["active_constraints"])
    println("  F:                  ", signif(rank_data[i]["flexibility_index"], 5))
end
