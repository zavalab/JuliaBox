
# Loading packages:
using JuMP 
using Distributions 
using Ipopt
using Gadfly
set_default_plot_size(20cm, 15cm)

# Generate random data: 
NS = 100
S = collect(1:NS); 
srand(0)
μ = 0; σ = 2; 
d = Normal(μ,σ)
R = rand(d,NS);

# Plotting data
p = plot(x = R, Geom.histogram(bincount = 30),Guide.XLabel("ξ"), Guide.YLabel("p(ξ)"))
display(p)

# ex1gauss.mod 
m = Model(solver=IpoptSolver(print_level=0))

@variable(m, x)            # decision variable
@variable(m, cost[1:NS])   # per scenario cost
@constraint(m, costeq[s in S], cost[s] == (x-R[s])^2) 

@objective(m, Min, (1/NS)*sum{cost[s], s in S})

solve(m)
solcost = getvalue(cost);
println("x=",getvalue(x))
println("mean=",getobjectivevalue(m))
println("meanR=",(1/NS)*sum(R))

# Plotting cost fistribution 
p = plot(x = solcost, Geom.histogram(bincount = 30),Guide.XLabel("f(x<sup>*</sup>(ξ),ξ)"), Guide.YLabel("p(f(⋅))"))
display(p)

m = Model(solver=IpoptSolver(print_level=0))

@variable(m, x)              # decision variable
@variable(m, cost[1:NS])     # per scenario cost
@variable(m, VaR)            # cvar auxiliary variable
@variable(m, phi[S] >= 0)    # cvar auxiliary variable

# set cvar probability level
alpha = 0.0001;               

@constraint(m, costeq[s in S], cost[s] == (x-R[s])^2) 
@constraint(m, cvar[s in S], cost[s]-VaR <= phi[s])

@objective(m, Min, VaR + (1/NS)*sum{(1/alpha)*phi[s], s in S})

solve(m)
solcost2 = getvalue(cost);
println("x=",getvalue(x))
println("CVaR=",getobjectivevalue(m))
println("mean=",(1/NS)*sum(solcost2))

# Plotting cost distribution
p = plot(x = solcost2, Geom.histogram(bincount = 30),Guide.XLabel("f(x<sup>*</sup>(ξ),ξ)"), Guide.YLabel("p(f(⋅))"))
display (p)

solcost
