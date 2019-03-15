using JuMP
using Ipopt

m = Model(solver=IpoptSolver(tol=1e-8))

@variable(m, 0 <= d1 <= 200);
@variable(m, 0 <= d2 <= 1000);
@variable(m, entropy)

@constraint(m, fix_demand, d1 + d2 == 800)
@NLconstraint(m, entropy == - (d1/(d1 + d2))*log(d1/(d1+d2)) - (d2/(d1 + d2))*log(d2/(d1+d2)))

@objective(m, Max, entropy)
solve(m)

println("Demand 1: ", getvalue(d1))
println("Demand 2: ", getvalue(d2))
println("Supply Node Price: ", getdual(fix_demand))
