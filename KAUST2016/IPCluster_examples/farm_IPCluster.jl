push!(LOAD_PATH, pwd())
include("ClusterIPM.jl")
using StochJuMP
using JuMP 
using Distributions 
using Ipopt


prcost = zeros(3)
prcost[1] = 150
prcost[2] = 230
prcost[3] = 260

pcost = zeros(3)
pcost[1] = 238
pcost[2] = 210
pcost[3] = 0

scost = zeros(3)
scost[1] = 170
scost[2] = 150
scost[3] = 36

demand = zeros(3)
demand[1] = 200
demand[2] = 240
demand[3] = 0;

NS = 10000;                 # number of scenarios
S = collect(1:NS)          # scenario set
P = collect(1:3)           # set of crops (1=wheat,2=corn,3=beets)

# assign random data
yield = zeros(length(S),3)
yield[S,1] = 2.5
yield[S,2] = 3.0
srand(123)
μ = 20; σ = 5
d = Normal(μ,σ)
yield[S,3] = rand(d,NS)


#=
m = Model(solver=IpoptSolver())
@defVar(m, x[S,P] >= 0)    # acres devoted to crops
@defVar(m, y[S,P] >= 0)    # crops purchase
@defVar(m, w[S,P] >= 0)    # crops sold;
@defVar(m, cost[s in S])
@addConstraint(m, varcost[s in S], cost[s] == sum{prcost[j]*x[s,j] + pcost[j]*y[s,j] - scost[j]*w[s,j], j in P}) 
@addConstraint(m, cap[s in S], sum{x[s,j], j in P} <= 500)
@addConstraint(m, bal[s in S,j in P], yield[s,j]*x[s,j]+y[s,j]-w[s,j] >= demand[j]) 
@addConstraint(m, sellb[s in S], w[s,3] <= 6000)
@addConstraint(m, buyb[s in S], y[s,3] <= 0)
@addConstraint(m, nonant[s in S,j in P], x[1,j] == x[s,j])
@setObjective(m, Min, (1/NS)*sum{cost[s], s in S})
solve(m)
=#


m = StochasticModel(NS)
@defVar(m, x[P] >= 0)    # acres devoted to crops
@defVar(m, s2 >= 0)
@addConstraint(m, cap, sum{x[j], j in P} + s2 == 500)
@setObjective(m, Min, sum{prcost[j]*x[j], j in P})
for i in 1:NS
    bl = StochasticBlock(m)
    @defVar(bl, y[P] >= 0)    # crops purchase
    @defVar(bl, w[P] >= 0)    # crops sold;
    @defVar(bl, s[P] >= 0)
    @addConstraint(bl, bal[j in P], yield[i,j]*x[j]+y[j]-w[j] - s[j] == demand[j])
    setUpper(w[3], 6000)
    setUpper(y[3], 0)
    @setObjective(bl, Min, 1.0/NS*sum{pcost[j]*y[j] - scost[j]*w[j], j in P})
end
CluIPM_solve(m)