using JuMP
using Ipopt

m = Model(solver=IpoptSolver(tol=1e-8))

beta = 2

@variable(m, 0 <= d1 <= 200);
@variable(m, 0 <= d2 <= 1000);
@variable(m, ge_0)

@constraint(m, fix_demand, d1 + d2 == 800)
@NLconstraint(m, ge_0 == -1/2*(log(d1) + log(*d2) - 2* log((d1 + d2)/2 )) )

@objective(m, Min, ge_0)
solve(m)

println("Demand 1: ", getvalue(d1))
println("Demand 2: ", getvalue(d2))
println("Supply Node Price: ", getdual(fix_demand))
