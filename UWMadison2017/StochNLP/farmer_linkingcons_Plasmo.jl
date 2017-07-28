# Example on using linking constraints in PLASMO
# Yankai Cao and Victor M. Zavala
# University of Wisconsin-Madison, 2016


using JuMP
using Distributions
using Ipopt
using Plasmo
MPI.Init()  # Initialize MPI

srand(123)
NS = 100;                   # number of scenarios
NP = 1000;
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

# create plasmo model
graph = GraphModel()
master = Model()
master_node = add_node(graph,master)

@variable(master, x[P] >= 0)    # acres devoted to crops
@variable(master, s2 >= 0)
@constraint(master, cap, (sum(x[j] for j in P) + s2) == 200)
@objective(master, Min, sum(prcost[j]*x[j] for j in P))

child_nodes = Array{Plasmo.NodeOrEdge}(NS)
for i in 1:NS
    bl = Model()
    child_node = add_node(graph,bl)
    child_nodes[i] = child_node
    @variable(bl, y[P] >= 0)    # crops purchase
    @variable(bl, 0<=w[j in P] <= sellub[j in P])    # crops sold
    @variable(bl, s[P] >= 0)
    @linkconstraint(graph, [j in P], yield[i,j]*master_node[:x][j]+y[j]-w[j] - s[j] == demand[j])

    @variable(bl, cost)
    @constraint(bl, cost ==sum(pcost[j]*y[j] - scost[j]*w[j] for j in P))
    @objective(bl, Min, 1.0/NS*cost)
end

# impose expected value constraint on cost
@linkconstraint(graph, sum(getnode(graph,i)[:cost] for i in 2:NS + 1) >= 50000)

#Solve with Ipopt
#graph.solver = IpoptSolver()
#solve(graph)


# solve with PIPS-NLP
#pipsnlp_solve(graph,master_node,child_nodes)

# access solution
#println(getvalue(getvariable(m, :x)))
#println(getvalue(getvariable(getNode(m,"s1"), :w)))
#println(getvalue(getvariable(getNode(m,"s1"), :cost)))
#println(getvalue(getvariable(getNode(m,"s2"), :cost)))

MPI.Finalize()
