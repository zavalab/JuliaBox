using JuMP
using Ipopt
#using Gurobi

m = Model(solver=IpoptSolver(tol=1e-8))

@variable(m, 0 <= d1 <= 200);
@variable(m, 0 <= d2 <= 1000);
alpha = 3
@variable(m, alpha_fair)

@constraint(m, fix_demand, d1 + d2 <= 800)
# @constraint(m, demand1, d1 == f1)
# @constraint(m, demand2, d2 == f2)
@NLconstraint(m, alpha_fair == 1/(1 - alpha)*(d1^(1 - alpha) + d2^(1 - alpha) ))
# @constraint(m, nash == alpha1*d1 + alpha2*d2)

@objective(m, Max, alpha_fair)
solve(m)

println("Demand 1: ", getvalue(d1))
println("Demand 2: ", getvalue(d2))
println("Alpha Value: ", alpha)
