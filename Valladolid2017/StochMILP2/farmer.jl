using JuMP 
using Ipopt 

# Model parameters

NS = 5;                  # number of scenarios
S = collect(1:NS);       # scenario set
P = collect(1:3);        # set of crops (1=wheat,2=corn,3=beets)

# Data 
prcost = zeros(3)          # production (planting) cost
prcost[1] = 150;
prcost[2] = 230;
prcost[3] = 260;

pcost = zeros(3)           # purchase cost
pcost[1] = 238;
pcost[2] = 210;
pcost[3] = 0;

scost = zeros(3)           # sales cost
scost[1] = 170;
scost[2] = 150;
scost[3] = 36;

demand = zeros(3)          # demand
demand[1] = 200;
demand[2] = 240;
demand[3] = 0;

# assign random data
yield = zeros(length(S),3) # yields
yield[S,1] = 2.5;
yield[S,2] = 3.0;
yield[1,3] = 10;
yield[2,3] = 15;
yield[3,3] = 20;
yield[4,3] = 25;
yield[5,3] = 30;

# Create Model
m = Model(solver=IpoptSolver(print_level=0))

@variable(m, x[P] >= 0)      # acres devoted to crops
@variable(m, y[S,P] >= 0)    # crops purchase
@variable(m, w[S,P] >= 0)    # crops sold;
@variable(m, cost[s in S])   # per scenario cost

@constraint(m, varcost[s in S], cost[s] == sum(prcost[j]*x[j] + pcost[j]*y[s,j] - scost[j]*w[s,j] for  j in P)) 
@constraint(m, cap, sum(x[j] for  j in P) <= 500)
@constraint(m, bal[s in S,j in P], yield[s,j]*x[j]+y[s,j]-w[s,j] >= demand[j]) 
@constraint(m, sellb[s in S], w[s,3] <= 6000)
@constraint(m, buyb[s in S], y[s,3] <= 0)

@objective(m, Min, (1/NS)*sum(cost[s] for s in S))

solve(m)

# Results
println(getvalue(x))
println("")
println(getvalue(cost))
println("")
println("obj ",getobjectivevalue(m))

