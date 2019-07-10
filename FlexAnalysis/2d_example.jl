using FlexJuMP, JuMP, Gurobi

# Define the distribution params
covar = [2 1; 1 3]
means = [4; 5]
dev = 3 * sqrt.(diag(covar))

# Setup the model for design A
ma = FlexibilityModel(solver = GurobiSolver(OutputFlag = 0))
@randomvariable(ma, θ[i = 1:2], mean = means[i])
@constraint(ma, θ[1] + θ[2] - 14 <= 0)
@constraint(ma, θ[1] - 2θ[2] - 2 <= 0)
@constraint(ma, -θ[1] <= 0)
@constraint(ma, -θ[2] <= 0)
setuncertaintyset(ma, :Ellipsoid, covar)

# Setup the model for design B
mb = FlexibilityModel(solver = GurobiSolver(OutputFlag = 0))
@randomvariable(mb, θ[i = 1:2], mean = means[i])
@constraint(mb, 0.75θ[1] + θ[2] - 14 <= 0)
@constraint(mb, 0.75θ[1] - 2θ[2] - 2 <= 0)
@constraint(mb, -θ[1] <= 0)
@constraint(mb, -θ[2] <= 0)
setuncertaintyset(mb, :Ellipsoid, covar)

# Get the ellipsoidal results
status = solve(ma, active_constr = true)
status = solve(mb, active_constr = true)
F_ellip_a = getflexibilityindex(ma)
F_ellip_b = getflexibilityindex(mb)
α_a = getconfidencelevel(ma)
α_b = getconfidencelevel(mb)

# Get the hyperbox results
setuncertaintyset(ma, :Hyperbox, [[dev]; [dev]])
setuncertaintyset(mb, :Hyperbox, [[dev]; [dev]])
status = solve(ma, active_constr = true)
status = solve(mb, active_constr = true)
F_box_a = getflexibilityindex(ma)
F_box_b = getflexibilityindex(mb)

# Get the samples SF indexes
SF_a = findstochasticflexibility(ma, num_pts = 100000)
SF_b = findstochasticflexibility(mb, num_pts = 100000)

# Print the results
println("****Design Compare Results****")
println("Design A")
println("  F_box:   ", signif(F_box_a, 3))
println("  F_ellip: ", signif(F_ellip_a, 3))
println("  α*:      ", signif(α_a * 100, 3), "%")
println("  SF-MC:   ", signif(SF_a * 100, 3), "%")
println("Design B")
println("  F_box:   ", signif(F_box_b, 3))
println("  F_ellip: ", signif(F_ellip_b, 3))
println("  α*:      ", signif(α_b * 100, 3), "%")
println("  SF-MC:   ", signif(SF_b * 100, 3), "%")
println("")

# Rank the constraints of design A and print
setuncertaintyset(ma, :Ellipsoid, covar)
rank_data = rankinequalities(ma, active_constr = true)
println("****Ranking Results****")
for i = 1:length(rank_data)
    println("Rank ", i)
    println("  Active Constraint:    ", rank_data[i]["active_constraints"])
    println("  F_ellip:              ", signif(rank_data[i]["flexibility_index"], 3))
end
