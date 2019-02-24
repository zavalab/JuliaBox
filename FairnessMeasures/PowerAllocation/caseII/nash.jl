using JuMP
using Ipopt
using Gurobi

m = Model(solver=IpoptSolver(tol=1e-8))

@variable(m, 0 <= p1 <= 1000);
@variable(m, 0 <= p2 <= 1000);
@variable(m, nash)

@constraint(m, fix_demand, p1 + p2 <= 800)
@NLconstraint(m, nash == log(p1) + log(p2))

@objective(m, Max, nash)
solve(m)

println("Demand 1: ", getvalue(p1))
println("Demand 2: ", getvalue(p2))
println("Supply Node Price: ", getdual(fix_demand))
