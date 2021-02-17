using JuMP, InfiniteOpt, Distributions, Gurobi, Random, PyPlot

# seed the uncertainty
Random.seed!(42)

# Set the problem parameters
A = [-1 -1 -1 0 0; 1 0 0 1 0; 0 1 0 -1 -1; 0 0 1 0 1]
Cg, Cξ = [1 0; 0. 0; 0 0; 0 1], [0 0; 1 0; 0 1.; 0 0]
c = [1, 10]
yg_lim = 10 * ones(size(Cg,2))
yb_lim = 4 * ones(size(A, 2))
αs = [0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.15, 0.17, 0.2, 0.5, 1]
M = 100
μ = [3., 5.]
Σ = [2 0.; 0. 2]

# Define the first model (traditional joint-chance)
model1 = InfiniteModel(Gurobi.Optimizer)
set_silent(model1)
@finite_parameter(model1, α1, 0.05)
@infinite_parameter(model1, ξ[1:2] ∈ MvNormal(μ, Σ), num_supports = 1000)
@infinite_variable(model1, yg[1:2](ξ) >= 0)
@infinite_variable(model1, yb[1:5](ξ))
@infinite_variable(model1, w(ξ), Bin)
@objective(model1, Min, expect(c' * yg, ξ))
@constraint(model1, A * yb .+ Cg * yg .- Cξ * ξ .== 0)
@constraint(model1, yg - yg_lim .<= w * M)
@constraint(model1, -yb - yb_lim .<= w * M)
@constraint(model1, yb - yb_lim .<= w * M)
chance1 = expect(w, ξ)
@constraint(model1, chance1 <= α1)

# Define the second model (using the alternative constraint logic)
model2 = InfiniteModel(Gurobi.Optimizer)
set_silent(model2)
@finite_parameter(model2, α2, 0.05)
@infinite_parameter(model2, ξ[1:2] ∈ MvNormal(μ, Σ), num_supports = 1000)
@infinite_variable(model2, yg[1:2](ξ) >= 0)
@infinite_variable(model2, yb[1:5](ξ))
@infinite_variable(model2, w[[:g, :b, :o]](ξ), Bin)
@objective(model2, Min, expect(c' * yg, ξ))
@constraint(model2, A * yb .+ Cg * yg .- Cξ * ξ .== 0)
@constraint(model2, yg - yg_lim .<= w[:g] * M)
@constraint(model2, -yb - yb_lim .<= w[:b] * M)
@constraint(model2, yb - yb_lim .<= w[:b] * M)
@constraint(model2, w[:o] >= w[:g] + w[:b] - 1)
chance2 = expect(w[:o], ξ)
@constraint(model2, chance2 <= α2)

# Setup storage containers
probs1 = zeros(length(αs))
objs1 = zeros(length(αs))
probs2 = zeros(length(αs))
objs2 = zeros(length(αs))

# Solve for the Pareto frontier
for (i, prob) in enumerate(αs)
    set_value(α1, prob)
    set_value(α2, prob)
    optimize!(model1)
    if primal_status(model1) == MOI.FEASIBLE_POINT
        objs1[i] = objective_value(model1)
        probs1[i] = value(chance1)
    else
        println("Model 1 infeasible with α = ", prob)
    end 
    optimize!(model2)
    if primal_status(model2) == MOI.FEASIBLE_POINT
        objs2[i] = objective_value(model2)
        probs2[i] = value(chance2)
    else
        println("Model 2 infeasible with α = ", prob)
    end 
end 

# Add extra solution to model 2
set_value(α2, 0)
optimize!(model2)
objs2 = [objective_value(model2); objs2]
probs2 =  [value(chance2); probs2]

# Make a plot
figure()
plot(objs1, 1 .- probs1, "C0--")
scatter(objs1, 1 .- probs1, label = "Full Intersection", color = "C0")
plot(objs2, 1 .- probs2, "C1--")
scatter(objs2, 1 .- probs2, label = "Partial Intersection", color = "C1")
legend(loc = "best")
xlabel("Generation Cost")
ylabel(L"Probability Level ($1-\alpha$)")
# savefig("sopf_pareto.png")
