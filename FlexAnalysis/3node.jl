using FlexJuMP, JuMP, Gurobi, Pavito, Ipopt

# Specify solution parameters
center = :analytic # :analytic or :feasible
β = 0 # -40, 0, or 50
set_type = :Ellipsoid # :Ellipsoid or :Hyperbox
positive_set = false

# Setup the uncertainty set parameters
covar = [80. β β; β 80. β; β β 120.]
box_dev = 3 * sqrt.(diag(covar))

# Specify the network details
line_cap1 = [35; 50]
gen_cap1 = 100
line_cap2 = [17; 30]
gen_cap2 = [25; 50; 25]
line_cap3 = [35; 50; 20]
gen_cap3 = 100

# Setup the model for Design 1
m1 = FlexibilityModel(solver = PavitoSolver(mip_solver = GurobiSolver(OutputFlag = 0),
                     cont_solver = IpoptSolver(print_level = 0), log_level = 0, mip_solver_drives = false))

# Define variables
@randomvariable(m1, d[i = 1:3], mean = 0)
@recoursevariable(m1, a[1:2])
@recoursevariable(m1, s)

# Set the line capacity constraints
@constraint(m1, [line = 1:2], -line_cap1[line] <= a[line])
@constraint(m1, [line = 1:2], a[line] <= line_cap1[line])

# Set the generator capacity constraints
@constraint(m1, 0.0 <= s)
@constraint(m1, s <= gen_cap1)

# Set the node balance constraints
@constraint(m1, a[1] - d[1] == 0)
@constraint(m1, -a[1] -a[2] + s - d[2] == 0)
@constraint(m1, a[2] - d[3] == 0)

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
@randomvariable(m2, d[i = 1:3], mean = new_mean[i])
@recoursevariable(m2, a[1:2])
@recoursevariable(m2, s[1:3])

# Set the line capacity constraints
@constraint(m2, [line = 1:2], -line_cap2[line] <= a[line])
@constraint(m2, [line = 1:2], a[line] <= line_cap2[line])

# Set the generator capacity constraints
@constraint(m2, [gen = 1:3], 0.0 <= s[gen])
@constraint(m2, [gen = 1:3], s[gen] <= gen_cap2[gen])

# Set the node balance constraints
@constraint(m2, a[1] - d[1] + s[1] == 0)
@constraint(m2, -a[1] -a[2] + s[2] - d[2] == 0)
@constraint(m2, a[2] - d[3] + s[3] == 0)

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
@randomvariable(m3, d[i = 1:3], mean = new_mean[i])
@recoursevariable(m3, a[1:3])
@recoursevariable(m3, s)

# Set the line capacity constraints
@constraint(m3, [line = 1:3], -line_cap3[line] <= a[line])
@constraint(m3, [line = 1:3], a[line] <= line_cap3[line])

# Set the generator capacity constraints
@constraint(m3, 0.0 <= s)
@constraint(m3, s <= gen_cap3)

# Set the node balance constraints
@constraint(m3, a[1] - a[3] - d[1] == 0)
@constraint(m3, -a[1] -a[2] + s - d[2] == 0)
@constraint(m3, a[2] + a[3] - d[3] == 0)

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
rank_data = rankinequalities(m1, max_ranks = 6)
println("****Ranking Results****")
for i = 1:length(rank_data)
    println("Rank ", i)
    println("  Active Constraints: ", rank_data[i]["active_constraints"])
    println("  F:                  ", signif(rank_data[i]["flexibility_index"], 5))
end
