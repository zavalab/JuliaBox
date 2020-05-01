using Plasmo
using Plots

# High Level ModelGraph
parentgraph = ModelGraph()
@node(parentgraph, n0)
@variable(n0, x0 >= 0)
@constraint(n0,x0 <= 2)
@objective(n0,Min,x0)

# Add a node to the model graph
graph2 = ModelGraph()
add_subgraph!(parentgraph,graph2)

@node(graph2,nodes[4:6])
@variable(nodes[4], x[1:2] >= 0)
@constraint(nodes[4],x[1] + x[2] >= 4)
@objective(nodes[4],Min,x[1])
@linkconstraint(parentgraph, x0 + x[2] >= 4)

@variable(nodes[5], x[1:2] >= 0)
@constraint(nodes[5],x[1] + x[2] >= 2)
@objective(nodes[5],Min,x[1])
@linkconstraint(parentgraph, x0 + x[2] >= 3)

@variable(nodes[6], x[1:2] >= 0)
@constraint(nodes[6],x[1] + x[2] >= 2)
@objective(nodes[6],Min,x[1])
@linkconstraint(parentgraph, x0 + x[2] >= 1)

#Optimize with GLPK
optimize!(graph1,GLPK.Optimizer)

n0.label = "n0"
nodes[4].label = "n4"
nodes[5].label = "n5"
nodes[6].label = "n6"

plt_graph = Plots.plot(parentgraph,node_labels = true,markersize = 60,labelsize = 30,linewidth = 4,subgraph_colors = true);
plt_matrix = Plots.spy(parentgraph,node_labels = true);
