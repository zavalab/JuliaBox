using JuMP
using Ipopt

m = Model(solver=IpoptSolver(tol=1e-8))

beta = 2

@variable(m, 0 <= p1 <= 1000);
@variable(m, 0 <= p2 <= 1000);
@variable(m, ge_1)

@constraint(m, fix_demand, p1 + p2 <= 800)
@NLconstraint(m, ge_1 == 1/(2*(p1 + p2)/2)*(p1*log(p1) + p2*log(p2) -  (p1 + p2)*log((p1 + p2)/2 )) )

@objective(m, Min, ge_1)
solve(m)

println("Demand 1: ", getvalue(p1))
println("Demand 2: ", getvalue(p2))
println("Supply Node Price: ", getdual(fix_demand))
