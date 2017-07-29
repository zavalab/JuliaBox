# Example on how to use JUMP+Ipopt and PLASMO+Ipopt
# Yankai Cao, Jordan Jalving, Victor M. Zavala
# University of Wisconsin-Madison, 2016

using JuMP
using Distributions
using Ipopt
using Plasmo

# set data
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

NS = 5;                    # number of scenarios
S = collect(1:NS)          # scenario set
P = collect(1:3)           # set of crops (1=wheat,2=corn,3=beets)

# assign random data
yield = zeros(length(S),3) # yields
yield[S,1] = 2.5;
yield[S,2] = 3.0;
yield[1,3] = 10;
yield[2,3] = 15;
yield[3,3] = 20;
yield[4,3] = 25;
yield[5,3] = 30

# construct problem with JuMP and solve using IPOPT
m = Model(solver=IpoptSolver())
@variable(m, x[P] >= 0)
@variable(m, y[S,P] >= 0)
@variable(m, w[S,P] >= 0)
@variable(m, cost[s in S])
@constraint(m, varcost[s in S], cost[s] == sum(prcost[j]*x[j] + pcost[j]*y[s,j] - scost[j]*w[s,j] for j in P))
@constraint(m, cap, sum(x[j] for j in P) <= 500)
@constraint(m, bal[s in S,j in P], yield[s,j]*x[j]+y[s,j]-w[s,j] >= demand[j])
@constraint(m, sellb[s in S], w[s,3] <= 6000)
@constraint(m, buyb[s in S], y[s,3] <= 0)
@objective(m, Min, (1/NS)*sum(cost[s] for s in S))
solve(m)
println(getvalue(x))
println(getvalue(w))

# construct problem with PLASMO and solve using Ipopt
graph = GraphModel()
master = Model()
master_node = add_node(graph,master)
# add variables, objective, and constraints to parent node (first-stage)
@variable(master, x[P] >= 0)
@constraint(master, cap, sum(x[j] for j in P)  <= 500)
@objective(master, Min, sum(prcost[j]*x[j] for j in P))

# add variables, objective, and constraints to child nodes (second-stage)
child_nodes = Array{Plasmo.NodeOrEdge}(NS)
for i in 1:NS
    bl = Model()
    # add children to graph
    child_node = add_node(graph,bl)
    child_nodes[i] = child_node
    @variable(bl, y[P] >= 0)
    @variable(bl, w[P] >= 0)
    setupperbound(w[3], 6000)
    setupperbound(y[3], 0)
    @objective(bl, Min, 1.0/NS*sum(pcost[j]*y[j] - scost[j]*w[j] for j in P))

    @linkconstraint(graph, [j in P], yield[i,j]*master_node[:x][j]+y[j]-w[j] >= demand[j])
end

# call Ipopt for solution
graph.solver = IpoptSolver()
solve(graph)

# access solution and display results
println(getvalue(x))
for i in 1:NS
    println(getvalue(getindex(child_nodes[i],:w))) 
end