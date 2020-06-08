include((@__DIR__)*"/ocp_model.jl")

using KaHyPar

#plot layout and matrix
plt_graph4 = plot(graph,layout_options = Dict(:tol => 0.1,:C => 2, :K => 4, :iterations => 500),linealpha = 0.2,markersize = 6);
plt_matrix4 = spy(graph);

#partition
hypergraph,hyper_map = gethypergraph(graph) #create hypergraph object based on graph
partition_vector = KaHyPar.partition(hypergraph,8,configuration = :connectivity,imbalance = 0.01)
partition = Partition(graph,partition_vector,hyper_map)
make_subgraphs!(graph,partition)

#plot partition
plt_graph5 = plot(graph,layout_options = Dict(:tol => 0.01, :iterations => 500),linealpha = 0.2,markersize = 6,subgraph_colors = true);
plt_matrix5 = spy(graph,subgraph_colors = true);


#combine subgraphs
combined_graph,reference_map = combine(graph,0)

#plot combined graph layout and matrix
plt_graph6 = plot(combined_graph,layout_options = Dict(:tol => 0.01,:iterations => 10),node_labels = true,markersize = 30,labelsize = 20,node_colors = true);
plt_matrix6 = spy(combined_graph,node_labels = true,node_colors = true);
