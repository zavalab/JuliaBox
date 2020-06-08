include((@__DIR__)*"/ocp_model.jl")

using KaHyPar
#partition
hypergraph,hyper_map = gethypergraph(graph) #create hypergraph object based on graph
partition_vector = KaHyPar.partition(hypergraph,8,configuration = :connectivity,imbalance = 0.01)
partition = Partition(graph,partition_vector,hyper_map)
make_subgraphs!(graph,partition)

#expand
subgraphs = getsubgraphs(graph)
expanded_subs = expand.(Ref(graph),subgraphs,Ref(2))

#plot expanded subgraphs
plt_graph7 = plot(graph,expanded_subs,layout_options = Dict(:tol => 0.01,:iterations => 1000),markersize = 6,linealpha = 0.2)
plt_matrix7 = spy(graph,expanded_subs)
