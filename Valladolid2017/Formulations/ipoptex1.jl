
using JuMP 
using Ipopt 

# Model 
m = Model(solver=IpoptSolver(tol = 1e-4, max_iter = 100,linear_solver ="mumps",mu_strategy="monotone"))

@variable(m, x >= 0)    
@variable(m, y >= 0)    

@NLconstraint(m, cons, x^2+y^2==1) 

@NLobjective(m, Min, (x-1)^2 + (x-y)^3)

# print model
println("\n")
print(m)
println("\n")

# solve model and get solution
solve(m)
println("x = ", getvalue(x), "\ny = ", getvalue(y))
