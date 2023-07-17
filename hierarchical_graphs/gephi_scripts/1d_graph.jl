using Plasmo
using DelimitedFiles, CSV, DataFrames
using Gurobi

############## LOAD IN DATA ###########################

line_data = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Lines.csv"))
bus_data = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Buses.csv"))
gen_data = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/All_Generators.csv"))
gen_r = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Generators_R.csv"))
gen_DA = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Generators_DAUC.csv"))
gen_ST = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Generators_STUC.csv"))

gen_data_conv = vcat(gen_DA, gen_ST)

gen_data_DA = vcat(gen_DA, gen_r)

D_DA = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/DA_bus_demand.csv"))
D_ST = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/ST_bus_demand_15min.csv"))
D_RT = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/RT_bus_demand_15min.csv"))

wind_DA = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Wind/DA.csv"))
wind_ST = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Wind/ST_15min.csv"))
wind_RT = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Wind/RT_15min.csv"))

solar_DA = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Solar/DA.csv"))
solar_ST = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Solar/ST_15min.csv"))
solar_RT = DataFrame(CSV.File((@__DIR__)*"/datasets/nrel118/Solar/RT_15min.csv"))

xi_DA = Dict()
xi_ST = Dict()
xi_RT = Dict()

wind_names = names(wind_DA)[2:end]
solar_names = names(solar_DA)[2:end]

for i in 1:length(wind_names)
    xi_DA[wind_names[i]] = wind_DA[:, wind_names[i]]
    xi_ST[wind_names[i]] = wind_ST[:, wind_names[i]]
    xi_RT[wind_names[i]] = wind_RT[:, wind_names[i]]
end

for i in 1:length(solar_names)
    xi_DA[solar_names[i]] = solar_DA[:, solar_names[i]]
    xi_ST[solar_names[i]] = solar_ST[:, solar_names[i]]
    xi_RT[solar_names[i]] = solar_RT[:, solar_names[i]]
end

B = line_data[:, "Susceptance"]

phi_c = 25
phi_o = 25
phi_u = 5000

########################## LOAD IN OTHER .jl FILES ##########################################

include((@__DIR__)*"/layer_construction.jl")
include((@__DIR__)*"/link_solutions.jl")
include((@__DIR__)*"/layer_construction_export_figures.jl")

########################## DEFINE FUNCTION FOR BUILDING 1 DAY PROBLEM #######################

bus_len = size(bus_data, 1)
line_len = size(line_data, 1)
node_num = bus_len + line_len

node_region = Vector{Int}(zeros(node_num))

for i in 1:bus_len
    if bus_data[i, "Region"] == "R1"
        node_region[i] = 1
    elseif bus_data[i, "Region"] == "R2"
        node_region[i] = 2
    elseif bus_data[i, "Region"] == "R3"
        node_region[i] = 3
    else
        println("$i has an incorrect region")
    end
end

for i in 1:line_len#(bus_len + 1):(bus_len + line_len)
    bus_from = line_data[i, "Bus from"]
    bus_region = bus_data[bus_data[:, "Bus Name"] .== bus_from, "Region"][1]

    if length(bus_region) != 2
        println("Iteration $i has an error")
    end
    if bus_region == "R1"
        node_region[i + bus_len] = 1
    elseif bus_region == "R2"
        node_region[i + bus_len] = 2
    elseif bus_region == "R3"
        node_region[i + bus_len] = 3
    else
        println("$i has an incorrect region")
    end
end

########################## ASSIGN NODES FOR DA PROBLEM ######################################

node_id = Vector{Int}([i for i in 1:(node_num * 25)])
region_DA = Vector{Int}(zeros(node_num * 25))
time = Vector{Int}(zeros(node_num * 25))
time_agg = Vector{Int}(zeros(node_num * 25))
for i in 1:25
    assignment_range = (1 + (i - 1) * node_num):(i * node_num)
    region_DA[assignment_range] .= node_region[:]
    time[assignment_range]   .= i
    time_agg[assignment_range] .= ceil(i/5)
end

DA_node_data = DataFrame(["id" => node_id, "Region" => region_DA, "TimePoint" => time, "TimeAgg" => time_agg])
#CSV.write((@__DIR__)*"/node_data_DA.csv", DA_node_data)

########################## ASSIGN NODES FOR ST PROBLEM ######################################

node_id_ST = Vector{Int}([i for i in 1:(node_num * 16)])
region_ST = Vector{Int}(zeros(node_num * 16))

for i in 1:16
    assignment_range = (1 + (i - 1) * node_num):(i * node_num)
    region_ST[assignment_range] .= node_region[:]
end

ST_node_data = DataFrame(["id" => node_id_ST, "Region" => region_ST])
#CSV.write((@__DIR__)*"/node_data_ST.csv", ST_node_data)

########################## ASSIGN NODES FOR HA PROBLEM ######################################

node_id_HA = Vector{Int}([i for i in 1:(node_num * 5)])
region_HA = Vector{Int}(zeros(node_num * 5))

for i in 1:5
    assignment_range = (1 + (i - 1) * node_num):(i * node_num)
    region_HA[assignment_range] .= node_region[:]
end

HA_node_data = DataFrame(["id" => node_id_HA, "Region" => region_HA])
#CSV.write((@__DIR__)*"/node_data_HA.csv", HA_node_data)

########################## ASSIGN NODES FOR 1D PROBLEM ######################################
node_num_1d = 25 * node_num + 16 * 8 * node_num + 96 * 5 * node_num
node_id_1d = Vector{Int}([i for i in 1:(node_num_1d)])
time_1d = Vector{Int}(zeros(node_num_1d))
region_1d = Vector{Int}(zeros(node_num_1d))
subproblem_1d = Vector{Int}(zeros(node_num_1d))

for i in 1:25
    a_range = (1 + (i - 1) * node_num):(i * node_num)
    region_1d[a_range] .= node_region[:]
    time_1d[a_range] .= ceil(i / 6)
    subproblem_1d[a_range] .= 1
    if i == 25
        time_1d[a_range] .= 4
    end
end

for i in 1:8
    for j in 1:16
        a_range = (1 + (25 + j - 1 + (i - 1) * 16) * node_num):((25 + j + (i - 1) * 16) * node_num)
        region_1d[a_range] .= node_region[:]
        time_1d[a_range] .= ceil(i / 2)
        subproblem_1d[a_range] .= 2
    end
end

for i in 1:96
    for j in 1:5
        a_range = (1 + (25 + 16 * 8 + j - 1 + (i - 1) * 5) * node_num):((25 + 16 * 8 + j + (i - 1) * 5) * node_num)
        region_1d[a_range] .= node_region[:]
        time_1d[a_range] .= ceil(i / 24)
        subproblem_1d[a_range] .= 3
    end
end


node_data_1d = DataFrame(["id" => node_id_1d, "Region" => region_1d, "TimePoint" => time_1d, "Subproblem" => subproblem_1d])
#CSV.write((@__DIR__)*"/node_data_1d.csv", node_data_1d)

########################## DEFINE NODES FOR AGGREGATION #######################
# aggregate the buses together: every 304 nodes will be a different number
aglevel1 = Vector{Int}(zeros(node_num_1d))
# aggregate the subproblmes together: every subgraph will be a different number
aglevel2 = Vector{Int}(zeros(node_num_1d))

for i in 1:Int(node_num_1d / node_num)
    t_range = (1 + (i - 1) * 304):(i * 304)
    aglevel1[t_range] .= i
end

t_range_DA = 1:(node_num * 25)
aglevel2[t_range_DA] .= 1
for i in 1:8
    t_range_ST = (1 + node_num * 25 + (i - 1) * 16 * node_num):(node_num * 25 + i * 16 * node_num)
    aglevel2[t_range_ST] .= i + 1
end

for i in 1:96
    t_range_HA = (1 + node_num * 25 + 8 * 16 * node_num + (i - 1) * 5 * node_num):(node_num * 25 + 8 * 16 * node_num + i * 5 * node_num)
    aglevel2[t_range_HA] .= i + 9
end


node_data_1d_agg = DataFrame(["id" => node_id_1d, "Region" => region_1d, "TimePoint" => time_1d, "Subproblem" => subproblem_1d, "aglevel1" => aglevel1, "aglevel2" => aglevel2])
#CSV.write((@__DIR__)*"/node_data_1d_agg.csv", node_data_1d_agg)


########################## DEFINE FUNCTION FOR BUILDING 1 DAY PROBLEM #######################
graph = OptiGraph()

graph_DAUC = OptiGraph()
day_num = 1

build_bus_over_time(graph_DAUC, 25, D_DA[:, 3:27], xi_DA, 1, gen_DA, gen_data_DA, gen_DA)

link_DA_over_time(graph_DAUC, gen_data_DA)
link_between_DA(graph_DAUC, [graph_DAUC], 1)

add_subgraph!(graph, graph_DAUC)

for i in 1:8
    local graph_STUC = OptiGraph()

    t_range = (3 + (i - 1) * 12 + (day_num - 1) * 96):(2 + 16 + (i - 1) * 12 + (day_num - 1) * 96)

    offset = (i - 1) * 12 + (day_num - 1) * 96
    build_bus_over_time(graph_STUC, 16, D_ST[:, t_range], xi_ST, 0.25, gen_ST, gen_data, gen_ST, offset = offset)
    link_ST_over_time(graph_STUC)
    add_subgraph!(graph, graph_STUC)

    link_DA_to_ST(graph, (1 + (i - 1) * 3):(1 + i * 3), i)
end

gen_com_HA = gen_data[gen_data[:, "Fuel"] .== "EMPTY", :]

for i in 1:96
    local graph_HAED = OptiGraph()
    build_bus_over_time(graph_HAED, 5, D_RT[:, (3 + (i - 1)):(2 + (i - 1) + 5)], xi_RT, 0.25, gen_com_HA, gen_data, gen_data_conv, offset = (i - 1), HA = true)
    add_subgraph!(graph, graph_HAED)
end

for i in 1:96
    link_over_HA(graph, i)
end

link_between_ST(graph)

using SparseArrays, LightGraphs
function get_edge_list(g)
    am = LightGraphs.adjacency_matrix(g)
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

cg, cm = Plasmo.clique_graph(graph)
g = cg.graph
edge_list = get_edge_list(g)

#CSV.write((@__DIR__)*"/edge_list.csv", edge_list)

cg_DA, cm = Plasmo.clique_graph(getsubgraphs(graph)[1])
g_DA = cg_DA.graph
edge_list_DA = get_edge_list(g_DA)
CSV.write((@__DIR__)*"/edge_list_DA.csv", edge_list_DA)


################################### 1 TIME POINT GRAPH #########################################

graph_118 = OptiGraph()
build_bus_over_time(graph_118, 1, D_DA[:, 3:3], xi_DA, 1, gen_DA, gen_data_DA, gen_DA)

using LightGraphs, SparseArrays
cg_118, cm_118 = Plasmo.clique_graph(graph_118)

g_118 = cg_118.graph
am_118 = LightGraphs.adjacency_matrix(g_118)

I_118, J_118, Z_118 = findnz(am_118)
edge_list_118 = Matrix{Int}(zeros(length(I_118), 2))

for i in 1:length(I_118)
    if J_118[i] <= I_118[i]
        edge_list_118[i, :] = [J_118[i], I_118[i]]
    end
end
edge_list_118 = edge_list_118[(edge_list_118[:, 1] .!= 0), :]

df_edges_118 = DataFrame(["Source" => edge_list_118[:, 1], "Target" => edge_list_118[:, 2]])
#CSV.write((@__DIR__)*"/edge_list_118_graph.csv", df_edges_118)

########################## ASSIGN NODES FOR 4hr PROBLEM ######################################
node_id = Vector{Int}([i for i in 1:node_num])
node_region
node_type = Vector{Any}(zeros(304))
node_type[1:size(bus_data, 1)] .= "Bus"
node_type[(1 + size(bus_data, 1)):(size(bus_data, 1) + size(line_data, 1))] .= "Line"

node_data_118 = DataFrame(["id" => node_id, "Region" => node_region, "Type" => node_type])
#CSV.write((@__DIR__)*"/node_data_118.csv", node_data_118)

################################### ST Graph #########################################

include((@__DIR__)*"/layer_construction_STHA.jl")
include((@__DIR__)*"/link_solutions_STHA.jl")

graph = OptiGraph()

graph_DAUC = OptiGraph()
day_num = 1

build_bus_over_time(graph_DAUC, 25, D_DA[:, 3:27], xi_DA, 1, gen_DA, gen_data_DA, gen_DA)

link_DA_over_time(graph_DAUC, gen_data_DA)
link_between_DA(graph_DAUC, [graph_DAUC], 1)

add_subgraph!(graph, graph_DAUC)

for i in 1:8
    local graph_STUC = OptiGraph()

    t_range = (3 + (i - 1) * 12 + (day_num - 1) * 96):(2 + 16 + (i - 1) * 12 + (day_num - 1) * 96)

    offset = (i - 1) * 12 + (day_num - 1) * 96
    build_bus_over_time(graph_STUC, 16, D_ST[:, t_range], xi_ST, 0.25, gen_ST, gen_data, gen_ST, offset = offset)
    link_ST_over_time(graph_STUC)
    add_subgraph!(graph, graph_STUC)

    link_DA_to_ST(graph, (1 + (i - 1) * 3):(1 + i * 3), i)
end

gen_com_HA = gen_data[gen_data[:, "Fuel"] .== "EMPTY", :]

for i in 1:96
    local graph_HAED = OptiGraph()
    build_bus_over_time(graph_HAED, 5, D_RT[:, (3 + (i - 1)):(2 + (i - 1) + 5)], xi_RT, 0.25, gen_com_HA, gen_data, gen_data_conv, offset = (i - 1), HA = true)
    add_subgraph!(graph, graph_HAED)
end

for i in 1:96
    link_over_HA(graph, i)
end

link_between_ST(graph)

using SparseArrays, LightGraphs
function get_edge_list(g)
    am = LightGraphs.adjacency_matrix(g)
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

ST_graph = getsubgraphs(graph)[2]

HA_graph = OptiGraph()
build_bus_over_time(HA_graph, 5, D_RT[:, 3:7], xi_RT, 0.25, gen_ST, gen_data, gen_ST, offset = 0)
link_ST_over_time(HA_graph)
cg_HA, cm = Plasmo.clique_graph(HA_graph)
g_HA = cg_HA.graph


cg_ST, cm = Plasmo.clique_graph(ST_graph)
g_ST = cg_ST.graph
cg_HA, cm = Plasmo.clique_graph(HA_graph)
g_HA = cg_HA.graph
edge_list_ST = get_edge_list(g_ST)
edge_list_HA = get_edge_list(g_HA)

#CSV.write((@__DIR__)*"/edge_list_ST.csv", edge_list_ST)
#CSV.write((@__DIR__)*"/edge_list_HA.csv", edge_list_HA)
