using GraphIO

lgraph,proj_map = clique_expansion(hypergraph)
# lgraph = graph.lightgraph

#Create edge list
savegraph("gas_network.csv",lgraph,"gas_network",GraphIO.EdgeList.EdgeListFormat())
f = open("gas_network_attributes.csv", "w")

write(f,"id,partition,component\n")

for pipe in pipelines
    for node in all_nodes(pipe)
        node.ext[:component] = 1
    end
end

for compressor in compressors
    for node in all_nodes(compressor)
        node.ext[:component] = 2
    end
end

for junction in junctions
    for node in all_nodes(junction)
        node.ext[:component] = 3
    end
end

for (id,node) in enumerate(all_nodes(gas_network))
    partition = node_vector[id]
    component = node.ext[:component]
    write(f,"$id,$partition,$component\n")
end
close(f)
