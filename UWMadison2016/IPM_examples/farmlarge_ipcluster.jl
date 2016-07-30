push!(LOAD_PATH, pwd())
using JuMP 
using Distributions 
using Ipopt
using Plasmo


srand(123)
NS = 1000;                    # number of scenarios
NP = 100;
S = collect(1:NS)           # scenario set
P = collect(1:NP)           # set of crops (1=wheat,2=corn,3=beets)

prcost = zeros(NP)
d = Uniform(100,200)
prcost = rand(d,NP)

pcost = zeros(NP)
d = Uniform(100,200)
pcost = rand(d,NP)

scost = zeros(NP)
scost = pcost - 50

demand = zeros(NP)
d = Uniform(100,300)
demand = rand(d,NP)/NP

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
#=
m = Model(solver=IpoptSolver())
@variable(m, x[S,P] >= 0)    # acres devoted to crops
@variable(m, y[S,P] >= 0)    # crops purchase
@variable(m, 0<=w[S,j in P] <= sellub[j in P])    # crops sold;
@variable(m, cost[s in S])
@constraint(m, varcost[s in S], cost[s] == sum{prcost[j]*x[s,j] + pcost[j]*y[s,j] - scost[j]*w[s,j], j in P}) 
@constraint(m, cap[s in S], sum{x[s,j], j in P} <= 200)
@constraint(m, bal[s in S,j in P], yield[s,j]*x[s,j]+y[s,j]-w[s,j] >= demand[j]) 
@constraint(m, nonant[s in S,j in P], x[1,j] == x[s,j])
@objective(m, Min, (1/NS)*sum{cost[s], s in S})
solve(m)
=#

m = NetModel()
@variable(m, x[P] >= 0)    # acres devoted to crops
@variable(m, s2 >= 0)
@constraint(m, cap, (sum{x[j], j in P} + s2) == 200)
@objective(m, Min, sum{prcost[j]*x[j], j in P})
for i in 1:NS
    bl = Model()
    @variable(bl, y[P] >= 0)    # crops purchase
    @variable(bl, 0<=w[j in P] <= sellub[j in P])    # crops sold
    @variable(bl, s[P] >= 0)
    @constraint(m, bal[j in P], yield[i,j]*x[j]+y[j]-w[j] - s[j] == demand[j])
    @objective(bl, Min, 1.0/NS*sum{pcost[j]*y[j] - scost[j]*w[j], j in P})
    @addNode(m, bl, "s$i")
end
ClusterIPM_solve(m)
println(getvalue(getvariable(m,:x)))
println(getvalue(getvariable(getNode(m,"s1"),:w)))



