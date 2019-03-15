using JuMP
using Ipopt

m = Model(solver=IpoptSolver(tol=1e-8))

alpha1 = 2
alpha2 = 2

@variable(m, 0 <= p1 <= 1000*alpha1);
@variable(m, 0 <= p2 <= 1000*alpha2);
@variable(m, entropy)

@constraint(m, fix_demand, p1/alpha1 + p2/alpha2 <= 800)
@NLconstraint(m, entropy == - (p1/(p1 + p2))*log(p1/(p1+p2)) - (p2/(p1 + p2))*log(p2/(p1+p2)))

@objective(m, Max, entropy)
solve(m)

println("Demand 1: ", getvalue(p1/alpha1))
println("Demand 2: ", getvalue(p2/alpha2))
println("Supply Node Price: ", getdual(fix_demand))
