using JuMP
using Gurobi

m = Model(solver=GurobiSolver())

@variable(m, 0 <= d1 <= 1000)
@variable(m, 0 <= d2 <= 1000)
@variable(m, sw)

@constraint(m, fix_demand, d1 + d2 <= 800)
@constraint(m, sw == d1 + d2)

@objective(m, Max, sw)
solve(m)

println("Demand 1: ", getvalue(d1))
println("Demand 2: ", getvalue(d2))
println("Supply Node Price: ", getdual(fix_demand))
