using FlexJuMP, JuMP, Gurobi, Pavito, Ipopt

# Specify solution parameters
center = :feasible
set_type = :Ellipsoid # :Ellipsoid or :Hyperbox
positive_set = true

# Setup the uncertainty set parameters
covar = eye(10) * 1200
box_dev = 3 * sqrt.(diag(covar))

# Specify the network details
line_cap1 = 100
gen_cap1 = [332; 140; 100; 100; 100]
line_cap2 = 100
gen_cap2 = [332; 140; 100; 100; 100]
line_cap3 = 200
gen_cap3 = [432; 340]

# Setup the model for Design 1
m1 = FlexibilityModel(solver = PavitoSolver(mip_solver = GurobiSolver(OutputFlag = 0),
                     cont_solver = IpoptSolver(print_level = 0), log_level = 0, mip_solver_drives = false))

# Define variables
@randomvariable(m1, d[i = 1:11, mean = 0) # Temperarily set the mean to 0
@recoursevariable(m1, a[1:20])
@recoursevariable(m1, g[1:5])

# Set the line capacity constraints
@constraint(m1, [line = 1:20], -line_cap1 <= a[line])
@constraint(m1, [line = 1:20], a[line] <= line_cap1)

# Set the generator capacity constraints
@constraint(m1, [gen = 1:5], 0.0 <= g[gen])
@constraint(m1, [gen = 1:5], g[gen] <= gen_cap1[gen])

# Set the node balance constraints
@constraint(m1, g[1] - a[1] - a[6] == 0)
@constraint(m1, a[1] + g[2] - sum(a[i] for i = [2; 4; 5]) - d[1] == 0)
@constraint(m1, g[3] + a[2] - a[3] - d[2] == 0)
@constraint(m1, sum(a[i] for i = [3; 4; 8]) - sum(a[i] for i = [7; 11]) - d[3] == 0)
@constraint(m1, sum(a[i] for i = [5; 6; 7; 12]) - d[4] == 0)
@constraint(m1, g[4] + sum(a[i] for i = [16; 18]) - sum(a[i] for i = [12; 19]) - d[5] == 0)
@constraint(m1, a[9] - sum(a[i] for i = [8; 10]) == 0)
@constraint(m1, g[5] - a[9] == 0)
@constraint(m1, sum(a[i] for i = [10; 11]) - sum(a[i] for i = [13; 14]) - d[6] == 0)
@constraint(m1, sum(a[i] for i = [13; 20]) - d[7] == 0)
@constraint(m1, a[19] - a[20] - d[8] == 0)
@constraint(m1, a[17] - a[18] - d[9] == 0)
@constraint(m1, a[15] - sum(a[i] for i = [16; 17]) - d[10] == 0)
@constraint(m1, a[14] - a[15] - d[11] == 0)

# Define the uncertainty set
if set_type == :Ellipsoid
    setuncertaintyset(m1, :Ellipsoid, covar, only_positive = positive_set)
else
    setuncertaintyset(m1, :Hyperbox, [[box_dev]; [box_dev]], only_positive = positive_set)
end

# Compute the center
new_mean = findcenteredmean(m1, center = center, update_mean = true, only_positive = positive_set)

# Setup the model for Design 2
m2 = FlexibilityModel(solver = PavitoSolver(mip_solver = GurobiSolver(OutputFlag = 0),
                     cont_solver = IpoptSolver(print_level = 0), log_level = 0, mip_solver_drives = false))

# Define variables
@randomvariable(m2, d[i = 1:11, mean = new_mean[i])
@recoursevariable(m2, a[1:20])
@recoursevariable(m2, g[1:5])

# Set the line capacity constraints
@constraint(m2, [line = 1:20], -line_cap2 <= a[line])
@constraint(m2, [line = 1:20], a[line] <= line_cap2)

# Set the generator capacity constraints
@constraint(m2, [gen = 1:5], 0.0 <= g[gen])
@constraint(m2, [gen = 1:5], g[gen] <= gen_cap2[gen])

# Set the node balance constraints
@constraint(m2, g[1] - a[1] - a[6] == 0)
@constraint(m2, a[1] + g[2] - sum(a[i] for i = [2; 4]) - d[1] == 0)
@constraint(m2, g[3] + a[2] - a[3] - d[2] == 0)
@constraint(m2, sum(a[i] for i = [3; 4; 8]) - sum(a[i] for i = [7; 11]) - d[3] == 0)
@constraint(m2, sum(a[i] for i = [6; 7; 12]) - d[4] == 0)
@constraint(m2, g[4] + sum(a[i] for i = [16; 18]) - sum(a[i] for i = [12; 19]) - d[5] == 0)
@constraint(m2, a[9] - sum(a[i] for i = [8; 10]) == 0)
@constraint(m2, g[5] - a[9] == 0)
@constraint(m2, sum(a[i] for i = [10; 11]) - sum(a[i] for i = [13]) - d[6] == 0)
@constraint(m2, sum(a[i] for i = [13; 20]) - d[7] == 0)
@constraint(m2, a[19] - a[20] - d[8] == 0)
@constraint(m2, a[17] - a[18] - d[9] == 0)
@constraint(m2, a[15] - sum(a[i] for i = [16; 17]) - d[10] == 0)
@constraint(m2, -a[15] - d[11] == 0)

# Define the uncertainty set
if set_type == :Ellipsoid
    setuncertaintyset(m2, :Ellipsoid, covar, only_positive = positive_set)
else
    setuncertaintyset(m2, :Hyperbox, [[box_dev]; [box_dev]], only_positive = positive_set)
end

# Setup the model for Design 3
m3 = FlexibilityModel(solver = PavitoSolver(mip_solver = GurobiSolver(OutputFlag = 0),
                     cont_solver = IpoptSolver(print_level = 0), log_level = 0, mip_solver_drives = false))

# Define variables
@randomvariable(m1, d[i = 1:11, mean = new_mean[i])
@recoursevariable(m1, a[1:20])
@recoursevariable(m1, g[1:2])

# Set the line capacity constraints
@constraint(m1, [line = 1:20], -line_cap3 <= a[line])
@constraint(m1, [line = 1:20], a[line] <= line_cap3)

# Set the generator capacity constraints
@constraint(m1, [gen = 1:2], 0.0 <= g[gen])
@constraint(m1, [gen = 1:2], g[gen] <= gen_cap3[gen])

# Set the node balance constraints
@constraint(m1, g[1] - a[1] - a[6] == 0)
@constraint(m1, a[1] + g[2] - sum(a[i] for i = [2; 4; 5]) - d[1] == 0)
@constraint(m1, a[2] - a[3] - d[2] == 0)
@constraint(m1, sum(a[i] for i = [3; 4; 8]) - sum(a[i] for i = [7; 11]) - d[3] == 0)
@constraint(m1, sum(a[i] for i = [5; 6; 7; 12]) - d[4] == 0)
@constraint(m1, sum(a[i] for i = [16; 18]) - sum(a[i] for i = [12; 19]) - d[5] == 0)
@constraint(m1, a[9] - sum(a[i] for i = [8; 10]) == 0)
@constraint(m1, -a[9] == 0)
@constraint(m1, sum(a[i] for i = [10; 11]) - sum(a[i] for i = [13; 14]) - d[6] == 0)
@constraint(m1, sum(a[i] for i = [13; 20]) - d[7] == 0)
@constraint(m1, a[19] - a[20] - d[8] == 0)
@constraint(m1, a[17] - a[18] - d[9] == 0)
@constraint(m1, a[15] - sum(a[i] for i = [16; 17]) - d[10] == 0)
@constraint(m1, a[14] - a[15] - d[11] == 0)

# Define the uncertainty set
if set_type == :Ellipsoid
    setuncertaintyset(m3, :Ellipsoid, covar, only_positive = positive_set)
else
    setuncertaintyset(m3, :Hyperbox, [[box_dev]; [box_dev]], only_positive = positive_set)
end

# Solve Design 1 and get the results
status1 = solve(m1)
F1 = getflexibilityindex(m1)
time1 = getsolutiontime(m1)
SF1 = findstochasticflexibility(m1, num_pts = 1000000, use_vulnerability_model = true, only_positive = positive_set)
if set_type == :Ellipsoid
    α1 = getconfidencelevel(m1)
end

# Solve Design 2 and get the results
status2 = solve(m2)
F2 = getflexibilityindex(m2)
time2 = getsolutiontime(m2)
SF2 = findstochasticflexibility(m2, num_pts = 1000000, use_vulnerability_model = true, only_positive = positive_set)
if set_type == :Ellipsoid
    α2 = getconfidencelevel(m2)
end

# Solve Design 3 and get the results
status3 = solve(m3)
F3 = getflexibilityindex(m3)
time3 = getsolutiontime(m3)
SF3 = findstochasticflexibility(m3, num_pts = 1000000, use_vulnerability_model = true, only_positive = positive_set)
if set_type == :Ellipsoid
    α3 = getconfidencelevel(m3)
end

# Print the results
println("****Comparison Results****")
println("Design 1")
println("  F index:  ", signif(F1, 4))
println("  Time:     ", signif(time1, 4))
println("  SF index: ", signif(SF1 * 100, 4), "%")
if set_type == :Ellipsoid
    println("  α*:       ", signif(α1 * 100, 4), "%")
end
println("Design 2")
println("  F index:  ", signif(F2, 4))
println("  Time:     ", signif(time2, 4))
println("  SF index: ", signif(SF2 * 100, 4), "%")
if set_type == :Ellipsoid
    println("  α*:       ", signif(α2 * 100, 4), "%")
end
println("Design 3")
println("  F index:  ", signif(F3, 4))
println("  Time:     ", signif(time3, 4))
println("  SF index: ", signif(SF3 * 100, 4), "%")
if set_type == :Ellipsoid
    println("  α*:       ", signif(α3 * 100, 4), "%")
end

# Rank the constraints of design A and print
rank_data = rankinequalities(m1, max_ranks = 30)
println("****Ranking Results****")
for i = 1:length(rank_data)
    println("Rank ", i)
    println("  Active Constraints: ", rank_data[i]["active_constraints"])
    println("  F:                  ", signif(rank_data[i]["flexibility_index"], 5))
end
