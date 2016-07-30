push!(LOAD_PATH, pwd())
using JuMP 
using Distributions 
using Ipopt
using Plasmo

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

NS = 100;                  # number of scenarios
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

# construct problem with JuMP and solve using IPOPT
m = Model(solver=IpoptSolver())
@variable(m, x[S,P] >= 0)    # acres devoted to crops
@variable(m, y[S,P] >= 0)    # crops purchase
@variable(m, w[S,P] >= 0)    # crops sold;
@variable(m, cost[s in S])
@constraint(m, varcost[s in S], cost[s] == sum{prcost[j]*x[s,j] + pcost[j]*y[s,j] - scost[j]*w[s,j], j in P})
@constraint(m, cap[s in S], sum{x[s,j], j in P} <= 500)
@constraint(m, bal[s in S,j in P], yield[s,j]*x[s,j]+y[s,j]-w[s,j] >= demand[j])
@constraint(m, sellb[s in S], w[s,3] <= 6000)
@constraint(m, buyb[s in S], y[s,3] <= 0)
@constraint(m, nonant[s in S,j in P], x[1,j] == x[s,j])
@objective(m, Min, (1/NS)*sum{cost[s], s in S})
solve(m)

# construct problem with PLASMO and solve using PIPSNLP
m = NetModel()
@variable(m, x[P] >= 0)    # acres devoted to crops
@constraint(m, cap, sum{x[j], j in P}  <= 500)
@objective(m, Min, sum{prcost[j]*x[j], j in P})
for i in 1:NS
    bl = Model()
    @variable(bl, y[P] >= 0)    # crops purchase
    @variable(bl, w[P] >= 0)    # crops sold;
    setupperbound(w[3], 6000)
    setupperbound(y[3], 0)
    @objective(bl, Min, 1.0/NS*sum{pcost[j]*y[j] - scost[j]*w[j], j in P})
    @addNode(m, bl, "s$i")
    @constraint(m, bal[j in P], yield[i,j]*x[j]+y[j]-w[j] >= demand[j])
end
ParPipsNlp_solve(m)
println(getvalue(getvariable(m, :x)))
println(getvalue(getvariable(getNode(m,"s1"), :w)))
println(getvalue(getvariable(getNode(m,"s2”), :w)))