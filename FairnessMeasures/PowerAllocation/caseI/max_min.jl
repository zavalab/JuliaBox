using JuMP
using Gurobi

m = Model(solver=GurobiSolver())

@variable(m, 0 <= d1 <= 200);
@variable(m, 0 <= d2 <= 1000);
@variable(m, phi)

@constraint(m, fix_demand, d1 + d2 <= 800)
@constraint(m, phi_cons1, phi <= d1)
@constraint(m, phi_cons2, phi <= d2)

@objective(m, Max, phi)
solve(m)

println("Demand 1: ", getvalue(d1))
println("Demand 2: ", getvalue(d2))
println("Supply Node Price: ", getdual(fix_demand))
