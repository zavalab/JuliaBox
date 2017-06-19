
# Loading packages:
using JuMP 
using Distributions 
using Ipopt
using PyPlot

# Generate random data: 
NS = 100
S = collect(1:NS); 
srand(0)
μ = 0; σ = 2; 
d = Normal(μ,σ)
R = rand(d,NS);

# Plotting data
plt[:hist](R, bins = 30);
grid("on")
xlabel(L"\xi")
ylabel(L"p(\xi)")

m = Model(solver=IpoptSolver(print_level=0))

@variable(m, x)            # decision variable
@variable(m, cost[1:NS])   # per scenario cost
@constraint(m, costeq[s in S], cost[s] == (x-R[s])^2 - R[s]*x) 
@objective(m, Min, (1/NS)*sum(cost[s] for s in S))

solve(m)
solcost = getvalue(cost);
println("x=",getvalue(x))
println("mean=",getobjectivevalue(m))

# Plotting cost fistribution 
plt[:hist](solcost,bins = 30)
grid("on")
axis([-10, 50, 0, 50])
xlabel(L"f(x^*(\xi),\xi)")
ylabel(L"p(f(\cdot))")

m = Model(solver=IpoptSolver(print_level=0))

@variable(m, x)              # decision variable
@variable(m, cost[1:NS])     # per scenario cost
@variable(m, VaR)            # cvar auxiliary variable
@variable(m, phi[S] >= 0)    # cvar auxiliary variable
alpha = 0.0001;              # cvar probability level

@constraint(m, costeq[s in S], cost[s] == (x-R[s])^2 - R[s]*x) 
@constraint(m, cvar[s in S], cost[s]-VaR <= phi[s])
@objective(m, Min, VaR + (1/NS)*sum((1/alpha)*phi[s] for s in S))

solve(m)
solcost2 = getvalue(cost);
println("x=",getvalue(x))
println("CVaR=",getobjectivevalue(m))
println("mean=",(1/NS)*sum(solcost2))

# Plotting cost distribution 
plt[:hist](solcost2,bins = 30)
grid("on")
axis([-10, 50, 0, 50])
xlabel(L"f(x^*(\xi),\xi)")
ylabel(L"p(f(\cdot))")


