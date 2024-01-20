using Revise
using PlasmoData, Graphs, MetaGraphs, Random, LinearAlgebra

dim3 = 1

Random.seed!(10)
random_array = rand(100, 100, dim3)

dg = matrix_to_graph(random_array)

dg_size = Base.summarysize(dg)

mg = MetaGraph(dg.g, 3.0)

for j in 1:100
    for i in 1:100
        set_prop!(mg, (j - 1) * 100 + i, :name, (i, j))
        for k in 1:dim3
            set_prop!(mg, (j - 1) * 100 + i, Symbol("weight$k"), random_array[i, j, k])
        end
    end
end

mg_size = Base.summarysize(mg)

println("The size of the DataGraph is ", dg_size, " bytes")
println("The size of the MetaGraph is ", mg_size, " bytes")
println()
println("The meta graph uses ", mg_size / dg_size, " times more storage for this problem than the DataGraph")
println()

dim3 = 10

Random.seed!(10)
random_array = rand(100, 100, dim3)

dg = matrix_to_graph(random_array)

dg_size = Base.summarysize(dg)

mg = MetaGraph(dg.g, 3.0)

for j in 1:100
    for i in 1:100
        set_prop!(mg, (j - 1) * 100 + i, :name, (i, j))
        for k in 1:dim3
            set_prop!(mg, (j - 1) * 100 + i, Symbol("weight$k"), random_array[i, j, k])
        end
    end
end

mg_size = Base.summarysize(mg)

println("The size of the DataGraph is ", dg_size, " bytes")
println("The size of the MetaGraph is ", mg_size, " bytes")
println()
println("The meta graph uses ", mg_size / dg_size, " times more storage for this problem than the DataGraph")
println()

dim3 = 25

Random.seed!(10)
random_array = rand(100, 100, dim3)

dg = matrix_to_graph(random_array)

dg_size = Base.summarysize(dg)

mg = MetaGraph(dg.g, 3.0)

for j in 1:100
    for i in 1:100
        set_prop!(mg, (j - 1) * 100 + i, :name, (i, j))
        for k in 1:dim3
            set_prop!(mg, (j - 1) * 100 + i, Symbol("weight$k"), random_array[i, j, k])
        end
    end
end

mg_size = Base.summarysize(mg)

println("The size of the DataGraph is ", dg_size, " bytes")
println("The size of the MetaGraph is ", mg_size, " bytes")
println()
println("The meta graph uses ", mg_size / dg_size, " times more storage for this problem than the DataGraph")
println()

dim3 = 100

Random.seed!(10)
random_array = rand(100, 100, dim3)

dg = matrix_to_graph(random_array)

dg_size = Base.summarysize(dg)

mg = MetaGraph(dg.g, 3.0)

for j in 1:100
    for i in 1:100
        set_prop!(mg, (j - 1) * 100 + i, :name, (i, j))
        for k in 1:dim3
            set_prop!(mg, (j - 1) * 100 + i, Symbol("weight$k"), random_array[i, j, k])
        end
    end
end

mg_size = Base.summarysize(mg)

println("The size of the DataGraph is ", dg_size, " bytes")
println("The size of the MetaGraph is ", mg_size, " bytes")
println()
println("The meta graph uses ", mg_size / dg_size, " times more storage for this problem than the DataGraph")
println()

println("Storage Breakdown: ")
println("nodes storage = ", Base.summarysize(dg.nodes), "; percent total =  ", Base.summarysize(dg.nodes) / dg_size, "%")
println()
println("node_map storage = ", Base.summarysize(dg.node_map), "; percent total =  ", Base.summarysize(dg.node_map) / dg_size, "%")
println()
println("node_data storage = ", Base.summarysize(dg.node_data), "; percent total =  ", Base.summarysize(dg.node_data) / dg_size, "%")
println()
println("edges storage = ", Base.summarysize(dg.edges), "; percent total =  ", Base.summarysize(dg.edges) / dg_size, "%")
println()
println("edge_map storage = ", Base.summarysize(dg.edge_map), "; percent total =  ", Base.summarysize(dg.edge_map) / dg_size, "%")
println()
println("edge_data storage = ", Base.summarysize(dg.edge_data), "; percent total =  ", Base.summarysize(dg.edge_data) / dg_size, "%")
println()
println("SimpleGraph storage = ", Base.summarysize(dg.g), "; percent total =  ", Base.summarysize(dg.g) / dg_size, "%")
println()
println("graph_data storage = ", Base.summarysize(dg.graph_data), "; percent total =  ", Base.summarysize(dg.graph_data) / dg_size, "%")
println()

# Data storage if we store edges using node names vs. node indices
edge_test = Vector{Tuple{Any, Any}}(undef, length(dg.edges))

for i in 1:length(dg.edges)
    edge_test[i] = (dg.nodes[dg.edges[i][1]], dg.nodes[dg.edges[i][2]])
end

println("If edges are stored with node names rather than node indices, it requires ", Base.summarysize(edge_test) / Base.summarysize(dg.edges), " times the memory")
