push!(LOAD_PATH, pwd())
include("ClusterIPM.jl")
using StochJuMP
using JuMP 
using Distributions 
using Ipopt


srand(123)

NS = 100;                    # number of scenarios
NP = 1000;
S = collect(1:NS)           # scenario set
P = collect(1:NP)           # set of crops (1=wheat,2=corn,3=beets)


prcost = zeros(NP)
d = Uniform(100,200)
prcost = rand(d,NP)

pcost = zeros(NP)
d = Uniform(100,200)
prcost = rand(d,NP)

scost = zeros(NP)
scost = prcost - 50

demand = zeros(NP)
d = Uniform(100,300)
prcost = rand(d,NP)

# assign random data
yield = zeros(length(S),NP)
d = Uniform(5,20)
for j in 1:(NP-1)
    yield[S,j] = rand(d,1)[1]
end
d = Uniform(10,30)
yield[S,NP] = rand(d,NS)

sellub = zeros(NP)
d = Uniform(2000,8000)
sellub[P] = rand(d,NP)

m = Model(solver=IpoptSolver())
@defVar(m, x[S,P] >= 0)    # acres devoted to crops
@defVar(m, y[S,P] >= 0)    # crops purchase
@defVar(m, 0<=w[S,j in P] <= sellub[j in P])    # crops sold;
@defVar(m, cost[s in S])
@addConstraint(m, varcost[s in S], cost[s] == sum{prcost[j]*x[s,j] + pcost[j]*y[s,j] - scost[j]*w[s,j], j in P}) 
@addConstraint(m, cap[s in S], sum{x[s,j], j in P} <= 200*NP)
@addConstraint(m, bal[s in S,j in P], yield[s,j]*x[s,j]+y[s,j]-w[s,j] >= demand[j]) 
@addConstraint(m, nonant[s in S,j in P], x[1,j] == x[s,j])
@setObjective(m, Min, (1/NS)*sum{cost[s], s in S})
solve(m)



m = StochasticModel(NS)
@defVar(m, x[P] >= 0)    # acres devoted to crops
@defVar(m, s2 >= 0)
@addConstraint(m, cap, sum{x[j], j in P} + s2 == 200*NP)
@setObjective(m, Min, sum{prcost[j]*x[j], j in P})
for i in 1:NS
    bl = StochasticBlock(m)
    @defVar(bl, y[P] >= 0)    # crops purchase
    @defVar(bl, 0<=w[j in P] <= sellub[j in P])    # crops sold;
    @defVar(bl, s[P] >= 0)
    @addConstraint(bl, bal[j in P], yield[i,j]*x[j]+y[j]-w[j] - s[j] == demand[j])
    @setObjective(bl, Min, 1.0/NS*sum{pcost[j]*y[j] - scost[j]*w[j], j in P})
end
CluIPM_solve(m)