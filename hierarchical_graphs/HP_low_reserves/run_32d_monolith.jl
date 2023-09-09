using Plasmo
using DelimitedFiles, CSV, DataFrames
using Gurobi

############## LOAD IN DATA ###########################

line_data = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Lines.csv"))
bus_data = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Buses.csv"))
gen_data = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/All_Generators.csv"))
gen_r = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Generators_R.csv"))
gen_DA = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Generators_DAUC.csv"))
gen_ST = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Generators_STUC.csv"))

gen_data_conv = vcat(gen_DA, gen_ST)

gen_data_DA = vcat(gen_DA, gen_r)

D_DA_no_reserves = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/DA_bus_demand_extrapolation.csv"))
D_ST_no_reserves = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/ST_bus_demand_15min_extrapolation.csv"))
D_RT_no_reserves = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/RT_bus_demand_15min_extrapolation.csv"))

wind_DA = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Wind/DA_extrapolation.csv"))
wind_ST = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Wind/ST_15min_extrapolation.csv"))
wind_RT = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Wind/RT_15min_extrapolation.csv"))

solar_DA = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Solar/DA_extrapolation.csv"))
solar_ST = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Solar/ST_15min_extrapolation.csv"))
solar_RT = DataFrame(CSV.File((@__DIR__)*"/../datasets/nrel118/Solar/RT_15min_extrapolation.csv"))

D_DA = copy(D_DA_no_reserves)
D_ST = copy(D_ST_no_reserves)
D_RT = copy(D_RT_no_reserves)

# Add the reserves to the demand values
D_DA[:, 3:end] .= D_DA_no_reserves[:, 3:end] .* 1.1
D_ST[:, 3:end] .= D_ST_no_reserves[:, 3:end] .* 1.1
D_RT[:, 3:end] .= D_RT_no_reserves[:, 3:end] .* 1.025

# Create dictionaries containing the renewable data
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


# Define costs
phi_c = 25 # curtailment cost
phi_o = 25 # overgeneration cost
phi_u = 5000 # unmet load cost

# Empty data frame for generators committed in the HA level
gen_com_HA = gen_data[gen_data[:, "Fuel"] .== "EMPTY", :]

########################## LOAD IN .jl FILES FOR BUILDING PROBLEM ##########################################

include((@__DIR__)*"/layer_construction.jl")
include((@__DIR__)*"/link_solutions.jl")

########################## DEFINE FUNCTION FOR BUILDING 1 DAY PROBLEM #######################

graph_set = []
# Define a function for building and solving the graph of the first day
function run_day1_monolith(graph_set)
    # Initialize full day graph
    graph = OptiGraph()

    # Initialize DAUC graph
    graph_DAUC = OptiGraph()
    day_num = 1

    # Add buses to the graph
    build_bus_over_time(graph_DAUC, 25, D_DA[:, 3:27], xi_DA, 1, gen_DA, gen_data_DA, gen_DA)
    # Add links across time for binary variable constraints and ramping constraints
    link_DA_over_time(graph_DAUC, gen_data_DA)
    # Add constraints
    link_between_DA(graph_DAUC, [graph_DAUC], 1)

    # Add DAUC subproblem to graph
    add_subgraph!(graph, graph_DAUC)

    # Build STUC subproblem graphs
    for i in 1:8
        local graph_STUC = OptiGraph()

        t_range = (3 + (i - 1) * 12 + (day_num - 1) * 96):(2 + 16 + (i - 1) * 12 + (day_num - 1) * 96)

        offset = (i - 1) * 12 + (day_num - 1) * 96
        # Add buses to graph
        build_bus_over_time(graph_STUC, 16, D_ST[:, t_range], xi_ST, 0.25, gen_ST, gen_data, gen_ST, offset = offset)
        # Link subproblem graph over time
        link_ST_over_time(graph_STUC)

        # Add STUC subgraph to full day graph
        add_subgraph!(graph, graph_STUC)

        # Link DAUC to STUC problem
        link_DA_to_ST(graph, (1 + (i - 1) * 3):(1 + i * 3), i)
    end

    # Define empty dataframe for HA layer
    gen_com_HA = gen_data[gen_data[:, "Fuel"] .== "EMPTY", :]

    # Build HAED subproblem graphs
    for i in 1:96
        local graph_HAED = OptiGraph()
        # Add buses to graph
        build_bus_over_time(graph_HAED, 5, D_RT[:, (3 + (i - 1)):(2 + (i - 1) + 5)], xi_RT, 0.25, gen_com_HA, gen_data, gen_data_conv, offset = (i - 1), HA = true)
        # Add HAED subgraph to full day graph
        add_subgraph!(graph, graph_HAED)
    end

    # Link between HA subproblem and other layers
    for i in 1:96
        link_over_HA(graph, i)
    end

    # Add up/down constraints and ramping constraints
    link_between_ST(graph)

    # Set optimizer and accompanying attributes
    set_optimizer(graph, Gurobi.Optimizer)
    set_optimizer_attribute(graph, "MIPGap", .05)
    set_optimizer_attribute(graph, "Threads", 4)
    set_optimizer_attribute(graph, "PreSparsify", 1)
    set_optimizer_attribute(graph, "LogFile", (@__DIR__)*"/gurobi_logs/logfile1.log")
    set_optimizer_attribute(graph, "NodefileStart", .0005)
    set_optimizer_attribute(graph, "TimeLimit", 36000)

    # Solve graph
    optimize!(graph)
    subgraph_set = getsubgraphs(graph)

    # Save solution
    for i in 1:length(subgraph_set)
        values = value.(all_variables(subgraph_set[i]))
        writedlm((@__DIR__)*"/results/32d_monolith/mipgap5_results$(i)_day1.csv", values, ',')
    end
    push!(graph_set, graph)
end

# Define a function for building and solving the graph of the second day
function run_day2_monolith(graph_set)
    ########################## Solve Day 2 ##################################
    # Initialize full day graph
    graph = OptiGraph()

    # Initialize DAUC graph
    graph_DAUC = OptiGraph()
    day_num = 2

    # Add buses to the graph
    build_bus_over_time(graph_DAUC, 25, D_DA[:, 27:51], xi_DA, 1, gen_DA, gen_data_DA, gen_DA, offset = 24)
    # Get set of DAUC subproblems
    graph_DAUC_set = []
    push!(graph_DAUC_set, getsubgraphs(graph_set[1])[1])
    push!(graph_DAUC_set, graph_DAUC)
    HA_subgraph_for_DA = getsubgraphs(graph_set[1])[105]

    # Add links across time for binary variable constraints and ramping constraints
    link_DA_over_time(graph_DAUC, gen_data_DA)
    # Add constraints
    link_between_DA(graph_DAUC, graph_DAUC_set, 2; HA_graph = HA_subgraph_for_DA)

    # Add DAUC subproblem to graph
    add_subgraph!(graph, graph_DAUC)

    # Build STUC subproblem graphs
    for i in 1:8
        local graph_STUC = OptiGraph()
        t_range = (3 + (i - 1) * 12 + (day_num - 1) * 96):(2 + 16 + (i - 1) * 12 + (day_num - 1) * 96)

        offset = (i - 1) * 12 + (day_num - 1) * 96
        # Add buses to graph
        build_bus_over_time(graph_STUC, 16, D_ST[:, t_range], xi_ST, 0.25, gen_ST, gen_data, gen_ST, offset = offset)
        # Link subproblem graph over time
        link_ST_over_time(graph_STUC)
        # Add STUC subgraph to full day graph
        add_subgraph!(graph, graph_STUC)

        # Link DAUC to STUC problem
        link_DA_to_ST(graph, (1 + (i - 1) * 3):(1 + i * 3), i)
    end

    # Define empty dataframe for HA layer
    gen_com_HA = gen_data[gen_data[:, "Fuel"] .== "EMPTY", :]

    # Build HAED subproblem graphs
    for i in 1:96
        local graph_HAED = OptiGraph()
        # Add buses to graph
        build_bus_over_time(graph_HAED, 5, D_RT[:, (3 + (i - 1) + 96):(2 + (i - 1) + 5 + 96)], xi_RT, 0.25, gen_com_HA, gen_data, gen_data_conv, offset = (i - 1 + 96), HA = true)
        # Add HAED subgraph to full day graph
        add_subgraph!(graph, graph_HAED)
    end

    # Link between HA subproblem and other layers
    for i in 1:96
        link_over_HA(graph, i)
    end

    # Add up/down constraints and ramping constraints
    link_between_ST(graph; day_num = 2, graph_last = graph_set[1])

    # Link previous day solutions to current day
    link_sol_over_HA_new_day(graph, graph_set[1])

    # Set optimizer and accompanying attributes
    set_optimizer(graph, Gurobi.Optimizer)
    set_optimizer_attribute(graph, "MIPGap", .05)
    set_optimizer_attribute(graph, "Threads", 4)
    set_optimizer_attribute(graph, "PreSparsify", 1)
    set_optimizer_attribute(graph, "LogFile", (@__DIR__)*"/gurobi_logs/logfile2.log")
    set_optimizer_attribute(graph, "NodefileStart", .0005)
    set_optimizer_attribute(graph, "TimeLimit", 36000)

    # Solve graph
    optimize!(graph)
    subgraph_set = getsubgraphs(graph)

    # Save solution
    for i in 1:length(subgraph_set)
        values = value.(all_variables(subgraph_set[i]))
        writedlm((@__DIR__)*"/results/32d_monolith/mipgap5_results$(i)_day2.csv", values, ',')
    end
    push!(graph_set, graph)
end

function run_dayi_monolith(graph_set, day)
    ########################## Solve Day i ##################################
    # Remove first entry of graph_set (storing many graphs in memory is expensive)
    if length(graph_set) >= 3
        popfirst!(graph_set)
    end

    graph_set_len = length(graph_set)

    # Initialize full day graph
    graph = OptiGraph()

    # Initialize DAUC graph
    graph_DAUC = OptiGraph()
    DA_range = (3 + 24 * (day - 1)):(3 + 24 * day)
    # Add buses to the graph
    build_bus_over_time(graph_DAUC, 25, D_DA[:, DA_range], xi_DA, 1, gen_DA, gen_data_DA, gen_DA, offset = (day - 1) * 24)
    # Get set of DAUC subproblems
    graph_DAUC_set = []
    for k in 1:length(graph_set)
        push!(graph_DAUC_set, getsubgraphs(graph_set[k])[1])
    end
    push!(graph_DAUC_set, graph_DAUC)
    HA_subgraph_for_DA = getsubgraphs(graph_set[graph_set_len])[105]

    # Add links across time for binary variable constraints and ramping constraints
    link_DA_over_time(graph_DAUC, gen_data_DA)
    # Add constraints
    link_between_DA(graph_DAUC, graph_DAUC_set, day; HA_graph = HA_subgraph_for_DA)

    # Add DAUC subproblem to graph
    add_subgraph!(graph, graph_DAUC)

    # Build STUC subproblem graphs
    for i in 1:8
        local graph_STUC = OptiGraph()
        t_range = (3 + (i - 1) * 12 + (day - 1) * 96):(2 + 16 + (i - 1) * 12 + (day - 1) * 96)

        offset = (i - 1) * 12 + (day - 1) * 96
        # Add buses to graph
        build_bus_over_time(graph_STUC, 16, D_ST[:, t_range], xi_ST, 0.25, gen_ST, gen_data, gen_ST, offset = offset)
        # Link subproblem graph over time
        link_ST_over_time(graph_STUC)
        # Add STUC subgraph to full day graph
        add_subgraph!(graph, graph_STUC)

        # Link DAUC to STUC problem
        link_DA_to_ST(graph, (1 + (i - 1) * 3):(1 + i * 3), i)
    end

    # Define empty dataframe for HA layer
    gen_com_HA = gen_data[gen_data[:, "Fuel"] .== "EMPTY", :]

    # Build HAED subproblem graphs
    for i in 1:96
        local graph_HAED = OptiGraph()
        t_range = (3 + (i - 1) + (day - 1) * 96):(2 + (i - 1) + 5 + (day - 1) * 96)
        offset = (i - 1 + (day - 1) * 96)
        # Add buses to graph
        build_bus_over_time(graph_HAED, 5, D_RT[:, t_range], xi_RT, 0.25, gen_com_HA, gen_data, gen_data_conv, offset = offset, HA = true)
        # Add HAED subgraph to full day graph
        add_subgraph!(graph, graph_HAED)
    end

    # Link between HA subproblem and other layers
    for i in 1:96
        link_over_HA(graph, i)
    end

    # Add up/down constraints and ramping constraints
    link_between_ST(graph; day_num = day, graph_last = graph_set[graph_set_len])

    # Link previous day solutions to current day
    link_sol_over_HA_new_day(graph, graph_set[graph_set_len])

    # Set optimizer and accompanying attributes
    set_optimizer(graph, Gurobi.Optimizer)
    set_optimizer_attribute(graph, "MIPGap", .05)
    set_optimizer_attribute(graph, "Threads", 8)
    set_optimizer_attribute(graph, "PreSparsify", 1)
    set_optimizer_attribute(graph, "LogFile", (@__DIR__)*"/gurobi_logs/logfile$(day).log")
    set_optimizer_attribute(graph, "NodefileStart", .0005)
    set_optimizer_attribute(graph, "TimeLimit", 36000)

    # Solve graph
    optimize!(graph)
    subgraph_set = getsubgraphs(graph)

    # Save solution
    for i in 1:length(subgraph_set)
        values = value.(all_variables(subgraph_set[i]))
        writedlm((@__DIR__)*"/results/32d_monolith/mipgap5_results$(i)_day$(day).csv", values, ',')
    end
    push!(graph_set, graph)
end


graph_set = []
run_times = []
# Solve day 1
d1_time = @elapsed run_day1_monolith(graph_set)
push!(run_times, d1_time)
# Solve day 2
d2_time = @elapsed run_day2_monolith(graph_set)
push!(run_times, d2_time)

# SOlve day 3-32
for i in 3:32
    di_time = @elapsed run_dayi_monolith(graph_set, i)
    push!(run_times, di_time)
end
