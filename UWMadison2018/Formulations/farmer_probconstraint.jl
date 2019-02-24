
using JuMP 
using Cbc
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
yield[5,3] = 30

# Model (with exact probabilistic constraint - mixed-integer formulation)

m = Model(solver=CbcSolver())

@variable(m, x[S,P] >= 0)     # acres devoted to crops
@variable(m, y[S,P] >= 0)     # crops purchase
@variable(m, w[S,P] >= 0)     # crops sold
@variable(m, cost[s in S])    # per scenario cost
@variable(m,probcost[s in S]) # per scenario cost for probabilistic constraint
@variable(m,wp[S], Bin)       # indicator variables for probabilistic constraint

alphap=NS # number of constraints allowed to be violated

@constraint(m, varcost[s in S], cost[s] == sum(prcost[j]*x[s,j] + pcost[j]*y[s,j] - scost[j]*w[s,j] for j in P)) 
@constraint(m, cap[s in S], sum(x[s,j] for j in P) <= 500)
@constraint(m, bal[s in S,j in 1:2], yield[s,j]*x[s,j]+y[s,j]-w[s,j] >= demand[j]) 
@constraint(m, probeq1[s in S], probcost[s] == yield[s,3]*x[s,3]+y[s,3]-w[s,3] - demand[3]) 
@constraint(m, probeq2[s in S], probcost[s] >= -wp[s]*6000) 
@constraint(m, probeq3[s in S], (1/NS)*sum(wp[s] for s in S) <= alphap/NS) 
@constraint(m, sellb[s in S], w[s,3] <= 6000)
@constraint(m, buyb[s in S], y[s,3] <= 0)
@constraint(m, nonant[s in S,j in P], x[1,j] == x[s,j])

@objective(m, Min, (1/NS)*sum(cost[s] for s in S))

solve(m)

# Results 
println(getvalue(cost))
println("")
println(getvalue(wp))
println("")
println(getvalue(probcost))
println("")
println("obj ",getobjectivevalue(m))

# Model (with cvar approximation)
mp = Model(solver=CbcSolver())

@variable(mp, x[S,P] >= 0)     # acres devoted to crops
@variable(mp, y[S,P] >= 0)     # crops purchase
@variable(mp, w[S,P] >= 0)     # crops sold;
@variable(mp, cost[s in S])    # per scenario cost
@variable(mp,probcost[s in S]) # per scenario cost for probabilistic constraint
@variable(mp,VaR)              # CVaR auxiliary variable
@variable(mp,CVaR)             # CVaR auxiliary variable
@variable(mp,phi[S]>=0)        # CVaR auxiliary variable

alpha=5/NS  # alpha largest constraints in CVaR

@constraint(mp, varcost[s in S], cost[s] == sum(prcost[j]*x[s,j] + pcost[j]*y[s,j] - scost[j]*w[s,j] for j in P)) 
@constraint(mp, cap[s in S], sum(x[s,j] for j in P) <= 500)
@constraint(mp, bal[s in S,j in 1:2], yield[s,j]*x[s,j]+y[s,j]-w[s,j] >= demand[j]) 
@constraint(mp, probeq1[s in S], probcost[s] == -(yield[s,3]*x[s,3]+y[s,3]-w[s,3] - demand[3])) 
@constraint(mp, probeq2[s in S], probcost[s] - VaR <= phi[s])
@constraint(mp, probeq3, CVaR == VaR + (1/NS)*(1/alpha)*sum(phi[s] for s in S)) 
@constraint(mp, probeq4, CVaR <= 0) 

@constraint(mp, sellb[s in S], w[s,3] <= 6000)
@constraint(mp, buyb[s in S], y[s,3] <= 0)
@constraint(mp, nonant[s in S,j in P], x[1,j] == x[s,j])

@objective(mp, Min, (1/NS)*sum(cost[s] for s in S))

solve(mp)

# Results 
println(getvalue(cost))
println("")
println(getvalue(probcost))
println("")
println(getvalue(CVaR))
println("")
println(getvalue(VaR))
println("")
println("obj ",getobjectivevalue(mp))


