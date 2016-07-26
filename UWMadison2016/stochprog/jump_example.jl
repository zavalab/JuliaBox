
# load necessary libraries
using JuMP
using Cbc

# declare model and solver
m = Model(solver=CbcSolver())

# add variables
@variable(m, 0 <= x <= 2 )
@variable(m, 0 <= y <= 30 )

# add objective
@objective(m, Max, 5x + 3*y )

# add constraints
@constraint(m, 1x + 5y <= 3.0 )

# print model in readable form
print(m)

# solve problem
status = solve(m)

# get optimal objective and solution
println("Objective value: ", getobjectivevalue(m))
println("x = ", getvalue(x))
println("y = ", getvalue(y))
