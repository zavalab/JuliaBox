using Plasmo
using Plots
using KaHyPar

N = 100  #number of intervals
T = 200  #number of time points
k = round(Int64,T/N) #time points per node
d = sin.(1:T)
d =  reshape(d, (k, div(length(d), k)))

graph = ModelGraph()
@node(graph,nodes[1:n_nodes])

#Node models
for (i,node) in enumerate(nodes)
    d_node = d[:,i]
    M = length(d_node)
    @variable(node, x[1:M+1])
    @variable(node, u[1:M])
    @constraint(node, dynamics[i in 1:M], x[i+1] == x[i] + u[i] + d[i])
    @objective(node, Min, 0.001*sum(x[i]^2 - 2*x[i]*d[i] for i in 1:M) + sum(u[i]^2 for i in 1:M))
end

#Link constraints
for i = 1:n_nodes - 1
    ni = getnode(graph,i)
    nj = getnode(graph,i+1)
    @linkconstraint(graph, ni[:x][end] == nj[:x][1])  #last state in partition i is first state in partition j
end

#First node gets initial condition
n1 = getnode(graph,1)
@constraint(n1,n1[:x][1] == 0)

plt_graph = plot(graph,layout_options = Dict(:tol => 0.04,:C => 0.5, :K => 2, :iterations => 500),markersize = 18,node_labels = true,labelsize = 12);
plt_matrix = spy(graph);

hypergraph,hyper_map = gethypergraph(graph) #create hypergraph object based on graph
partition_vector = KaHyPar.partition(hypergraph,8,configuration = :edge_cut)
partition = Partition(hypergraph,partition_vector,hyper_map)
ModelGraphs.make_subgraphs!(graph,partition,hyper_map)

plt_graph2 = plot(graph,layout_options = Dict(:tol => 0.01,:C => 2, :K => 4, :iterations => 500),markersize = 18,node_labels = true,labelsize = 12,subgraph_colors = true)
plt_matrix2 = spy(graph)

new_graph,reference = combine(graph,0)

plt_graph3 = plot(new_graph,layout_options = Dict(:tol => 0.01,:C => 2, :K => 4, :iterations => 500),node_labels = true,markersize = 30,labelsize = 20)
plt_matrix3 = spy(new_graph,node_labels = true)
