using Graphs, DelimitedFiles, CSV, DataFrames, SparseArrays

function get_edge_list(g)
    am = Graphs.adjacency_matrix(g)
    I_g, J_g, z_g = findnz(am)
    edge_list = Matrix{Int}(zeros(length(I_g), 2))
    for i in 1:length(I_g)
        if J_g[i] <= I_g[i]
            edge_list[i, :] .= [J_g[i], I_g[i]]
        end
    end
    edge_list = edge_list[(edge_list[:, 1] .!= 0), :]
    df = DataFrame(["Source" => edge_list[:, 1], "Target" => edge_list[:, 2]])
    return df
end

edge_list, header = readdlm((@__DIR__)*"/edge_list.csv", ',', header = true)
edge_list = Matrix{Int}(edge_list)
num_nodes = maximum(edge_list)

g = Graphs.SimpleGraph(Int(num_nodes))

for i in 1:size(edge_list, 1)
    Graphs.add_edge!(g, edge_list[i, 1], edge_list[i, 2])
    if i%100000 == 0
        println("done with $i")
    end
end

node_data = DataFrame(CSV.File((@__DIR__)*"/node_data_1d_agg.csv"))

g_agg1 = Graphs.SimpleGraph(633)
g_agg2 = Graphs.SimpleGraph(105)

for i in 1:size(edge_list, 1)
    node1_agg1 = node_data[edge_list[i, 1], "aglevel1"]
    node1_agg2 = node_data[edge_list[i, 1], "aglevel2"]
    node2_agg1 = node_data[edge_list[i, 2], "aglevel1"]
    node2_agg2 = node_data[edge_list[i, 2], "aglevel2"]
    if node1_agg1 != node2_agg1
        Graphs.add_edge!(g_agg1, node1_agg1, node2_agg1)
    end
    if node1_agg2 != node2_agg2
        Graphs.add_edge!(g_agg2, node1_agg2, node2_agg2)
    end
end

df_agg1 = get_edge_list(g_agg1)
df_agg2 = get_edge_list(g_agg2)

agg1_node_id = Vector{Int}([i for i in 1:633])
agg2_node_id = Vector{Int}([i for i in 1:105])
agg1_time = Vector{Int}(zeros(633))
agg1_subgraph = Vector{Int}(zeros(633))
agg2_subgraph = Vector{Int}(zeros(105))

agg1_time[1:6] .= 1
agg1_time[7:12] .= 2
agg1_time[13:18] .= 3
agg1_time[19:25] .= 4

for i in 1:8
    index = Int(ceil(i/2))
    t_range = (1 + 25 + (i - 1) * 16):(25 + i * 16)
    agg1_time[t_range] .= index
end

for i in 1:96
    index = Int(ceil(i / 24))
    t_range = (1 + 25 + 8 * 16 + (i - 1) * 5):(25 + 8 * 16 + i * 5)
    agg1_time[t_range] .= index
end

agg1_subgraph[1:25] .= 1
agg1_subgraph[26:(25 + 8 * 16)] .= 2
agg1_subgraph[(1 + 25 + 8 * 16):633] .= 3

agg2_subgraph[1] = 1
agg2_subgraph[2:9] .= 2
agg2_subgraph[10:105] .= 3

node_data_agg1 = DataFrame(["id" => agg1_node_id, "TimePoint" => agg1_time, "Subproblem" => agg1_subgraph])
node_data_agg2 = DataFrame(["id" => agg2_node_id, "Subproblem" => agg2_subgraph])
#CSV.write((@__DIR__)*"/node_data_agg1.csv", node_data_agg1)
#CSV.write((@__DIR__)*"/node_data_agg2.csv", node_data_agg2)
#CSV.write((@__DIR__)*"/edge_list_agg1.csv", df_agg1)
#CSV.write((@__DIR__)*"/edge_list_agg2.csv", df_agg2)

DA_offset = 25 * 304

node_data[:, "miniproblem"] .= 0
node_data[:, "minitimepoints"] .= 0

for i in 1:4
    node_vals = (1 + (i - 1) * 304):(i * 304)
    n_range = (1 + DA_offset + (i - 1) * 304):(DA_offset + i * 304)
    node_data[n_range, "miniproblem"] .= node_vals
    node_data[n_range, "minitimepoints"] .= 1
end

ST_offset = 16 * 8 * 304

for i in 1:4
    HA_offset = (i - 1) * 5 * 304
    for j in 1:5
        n_range = (1 + DA_offset + ST_offset + HA_offset + (j - 1) * 304):(DA_offset + ST_offset + HA_offset + j * 304)
        node_vals = (1 + (4 + (i - 1) * 5 + j - 1) * 304):((4 + (i - 1) * 5 + j) * 304)
        node_data[n_range, "miniproblem"] .= node_vals
        node_data[n_range, "minitimepoints"] .= 1 + i
    end
end

mini_g = Graph(24 * 304)

for i in 1:size(edge_list, 1)
    node1_mini = node_data[edge_list[i, 1], "miniproblem"]
    node2_mini = node_data[edge_list[i, 2], "miniproblem"]
    if (node1_mini != 0) && (node2_mini != 0)
        Graphs.add_edge!(mini_g, node1_mini, node2_mini)
    end
end

mini_edge_list = get_edge_list(mini_g)

mini_node_ids = [i for i in 1:7296]
mini_timepoints = node_data[node_data[:, "miniproblem"] .!= 0, "minitimepoints"]
mini_subproblems = node_data[node_data[:, "miniproblem"] .!= 0, "Subproblem"]

mini_node_data = DataFrame(["id" => mini_node_ids, "TimePoint" => mini_timepoints, "Subproblem" => mini_subproblems])

#CSV.write((@__DIR__)*"/miniproblem_node_data.csv", mini_node_data)
#CSV.write((@__DIR__)*"/miniproblem_edge_list.csv", mini_edge_list)

DA_offset = 25 * 304

node_data[:, "miniproblem2"] .= 0
node_data[:, "minitimepoints2"] .= 0
node_data[:, "nodesize2"] .= 0

for i in 1:2
    node_vals = (1 + (i - 1) * 304):(i * 304)
    n_range = (1 + DA_offset + (i - 1) * 304):(DA_offset + i * 304)
    node_data[n_range, "miniproblem2"] .= node_vals
    node_data[n_range, "minitimepoints2"] .= 1
    node_data[n_range, "nodesize2"] .= 100
end

ST_offset = 16 * 8 * 304

for i in 1:2
    HA_offset = (i - 1) * 5 * 304
    for j in 1:5
        n_range = (1 + DA_offset + ST_offset + HA_offset + (j - 1) * 304):(DA_offset + ST_offset + HA_offset + j * 304)
        node_vals = (1 + (2 + (i - 1) * 5 + j - 1) * 304):((2 + (i - 1) * 5 + j) * 304)
        node_data[n_range, "miniproblem2"] .= node_vals
        node_data[n_range, "minitimepoints2"] .= 1 + i
        node_data[n_range, "nodesize2"] .= 10
    end
end

mini_g2 = Graph(12 * 304)

for i in 1:size(edge_list, 1)
    node1_mini = node_data[edge_list[i, 1], "miniproblem2"]
    node2_mini = node_data[edge_list[i, 2], "miniproblem2"]
    if (node1_mini != 0) && (node2_mini != 0)
        Graphs.add_edge!(mini_g2, node1_mini, node2_mini)
    end
end

mini_edge_list2 = get_edge_list(mini_g2)

mini_node_ids2 = [i for i in 1:3648]
mini_timepoints2 = node_data[node_data[:, "miniproblem2"] .!= 0, "minitimepoints2"]
mini_subproblems2 = node_data[node_data[:, "miniproblem2"] .!= 0, "Subproblem"]
mini_node_size2 = node_data[node_data[:, "miniproblem2"] .!= 0, "nodesize2"]

mini_node_data2 = DataFrame(["id" => mini_node_ids2, "TimePoint" => mini_timepoints2, "Subproblem" => mini_subproblems2, "NodeSize" => mini_node_size2])

#CSV.write((@__DIR__)*"/miniproblem_node_data2.csv", mini_node_data2)
#CSV.write((@__DIR__)*"/miniproblem_edge_list2.csv", mini_edge_list2)



DA_offset = 25 * 304

node_data[:, "miniproblem3"] .= 0
node_data[:, "minitimepoints3"] .= 0
node_data[:, "nodesize3"] .= 0

for i in 1:5
    node_vals = (1 + (i - 1) * 304):(i * 304)
    n_range = (1 + DA_offset + (i - 1) * 304):(DA_offset + i * 304)
    node_data[n_range, "miniproblem3"] .= node_vals
    node_data[n_range, "minitimepoints3"] .= 1
    node_data[n_range, "nodesize3"] .= 100
end

ST_offset = 16 * 8 * 304

for i in 1:1
    HA_offset = (i - 1) * 5 * 304
    for j in 1:5
        n_range = (1 + DA_offset + ST_offset + HA_offset + (j - 1) * 304):(DA_offset + ST_offset + HA_offset + j * 304)
        node_vals = (1 + (5 + (i - 1) * 5 + j - 1) * 304):((5 + (i - 1) * 5 + j) * 304)
        node_data[n_range, "miniproblem3"] .= node_vals
        node_data[n_range, "minitimepoints3"] .= 1 + i
        node_data[n_range, "nodesize3"] .= 10
    end
end

mini_g3 = Graph(10 * 304)

for i in 1:size(edge_list, 1)
    node1_mini = node_data[edge_list[i, 1], "miniproblem3"]
    node2_mini = node_data[edge_list[i, 2], "miniproblem3"]
    if (node1_mini != 0) && (node2_mini != 0)
        Graphs.add_edge!(mini_g3, node1_mini, node2_mini)
    end
end

mini_edge_list3 = get_edge_list(mini_g3)

mini_node_ids3 = [i for i in 1:3040]
mini_timepoints3 = node_data[node_data[:, "miniproblem3"] .!= 0, "minitimepoints3"]
mini_subproblems3 = node_data[node_data[:, "miniproblem3"] .!= 0, "Subproblem"]
mini_node_size3 = node_data[node_data[:, "miniproblem3"] .!= 0, "nodesize3"]

mini_node_data3 = DataFrame(["id" => mini_node_ids3, "TimePoint" => mini_timepoints3, "Subproblem" => mini_subproblems3, "NodeSize" => mini_node_size3])

mini_edge_list3[:, "edge_color"] .= 0
for i in 1:size(mini_edge_list3, 1)
    node1 = mini_edge_list3[i, "Source"]
    node2 = mini_edge_list3[i, "Target"]
    sub1 = mini_node_data3[node1, "Subproblem"]
    sub2 = mini_node_data3[node2, "Subproblem"]

    if sub1 != sub2
        mini_edge_list3[i, "edge_color"] = 1
    else
        mini_edge_list3[i, "edge_color"] = sub1
    end
end

#CSV.write((@__DIR__)*"/miniproblem_node_data3.csv", mini_node_data3)
#CSV.write((@__DIR__)*"/miniproblem_edge_list3.csv", mini_edge_list3)
