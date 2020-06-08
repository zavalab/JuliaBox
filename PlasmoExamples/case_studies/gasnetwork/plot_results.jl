using Plots
pyplot()

#Plot the modelgraph representation of the problem
plt_graph1 = Plots.plot(gas_network,node_labels = false,markersize = 6,subgraph_colors = true,
linewidth = 1,linealpha = 0.3,layout_options = Dict(:tol => 0.01,:iterations => 1000));
savefig(plt_graph1,"13_pipeline_modelgraph.pdf")

#Plot the compressor policy
rank_zero = manager.mpi2j[0] #julia process representing rank 0
solution = fetch(@spawnat(rank_zero, pipsgraph))

#TODO: copy solution to actual modelnodes in a sensical way (in the solver).

#TODO: copy solution from combined graph to original network using combine_map
