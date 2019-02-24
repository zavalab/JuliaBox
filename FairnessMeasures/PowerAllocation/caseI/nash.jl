using JuMP
using Ipopt
using Gurobi

m = Model(solver=IpoptSolver(tol=1e-8))

@variable(m, 0 <= d1 <= 200);
@variable(m, 0 <= d2 <= 1000);
# @variable(m,  f1 >= 0)
# @variable(m,  f2 >= 0)
@variable(m, nash)

@constraint(m, fix_demand, d1 + d2 <= 800)
@NLconstraint(m, nash == log(d1) + log(d2))

@objective(m, Max, nash)
solve(m)

println("Demand 1: ", getvalue(d1))
println("Demand 2: ", getvalue(d2))
println("Supply Node Price: ", getdual(fix_demand))
