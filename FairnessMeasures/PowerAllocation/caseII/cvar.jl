using JuMP
using CPLEX

m = Model(solver=CplexSolver())


@variable(m, 0 <= d1 <= 1000);
@variable(m, 0 <= d2 <= 1000);
@variable(m, cvar)
@variable(m, risk)
@variable(m, phi1 >= 0)
@variable(m, phi2 >= 0)
@variable(m, disutility1)
@variable(m, disutility2)

alpha = 0.99
@constraint(m, disutility1 == 0 - d1)
@constraint(m, disutility2 == 0 - d2)

@constraint(m, fix_demand, d1 + d2 <= 800)
@constraint(m, phi1 >= disutility1 - risk)
@constraint(m, phi2 >= disutility2 - risk)
@constraint(m, cvar == risk + (1/2)*( phi1/(1-alpha) + phi2/(1-alpha)))

@objective(m, Min, cvar)
solve(m)

println("Demand 1: ", getvalue(d1))
println("Demand 2: ", getvalue(d2))
println("Alpha Value: ", alpha)
println("Disutility 1: ", getvalue(disutility1))
println("Disutility 2: ", getvalue(disutility2))
