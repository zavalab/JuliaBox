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

D_DA[:, 3:end] .= D_DA_no_reserves[:, 3:end] .* 1.05
D_ST[:, 3:end] .= D_ST_no_reserves[:, 3:end] .* 1.05
D_RT[:, 3:end] .= D_RT_no_reserves[:, 3:end] .* 1.0125

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

gen_com_HA = gen_data[gen_data[:, "Fuel"] .== "EMPTY", :]

########################## LOAD IN OTHER .jl FILES ##########################################

include((@__DIR__)*"/layer_construction.jl")
include((@__DIR__)*"/link_solutions.jl")

########################## DEFINE FUNCTION FOR BUILDING 1 DAY PROBLEM #######################

function run_day1_serial(graph_set)
    graph_DAUC = OptiGraph()

    build_bus_over_time(graph_DAUC, 25, D_DA[:, 3:27], xi_DA, 1, gen_DA, gen_data_DA, gen_DA)
    link_DA_over_time(graph_DAUC, gen_data_DA)
    link_between_DA(graph_DAUC, [graph_DAUC], 1)

    set_optimizer(graph_DAUC, Gurobi.Optimizer)

    set_optimizer_attribute(graph_DAUC, "MIPGap", .001)

    println("Optimizing DA Problem")
    optimize!(graph_DAUC)


    graph = OptiGraph()

    add_subgraph!(graph, graph_DAUC)
    set_optimizer(graph, Gurobi.Optimizer)

    DAUC_sols = value.(all_variables(graph_DAUC))
    writedlm((@__DIR__)*"/results/32d_serial/DAUC_results_day1.csv", DAUC_sols, ',')

    for i in 1:8
        local graph_STUC = OptiGraph()
        build_bus_over_time(graph_STUC, 16, D_ST[:, (3 + (i - 1) * 12):(2 + (i - 1) * 12 + 16)], xi_ST, 0.25, gen_ST, gen_data, gen_ST, offset = (i - 1) * 12)
        link_ST_over_time(graph_STUC)
        link_DA_sol_to_ST(graph, graph_STUC, (1 + (i - 1) * 3):(4 + (i - 1) * 3))
        set_optimizer(graph_STUC, Gurobi.Optimizer)
        set_optimizer_attribute(graph_STUC, "MIPGap", .005)
        add_subgraph!(graph, graph_STUC)
    end

    for i in 1:96
        local graph_HAED = OptiGraph()
        build_bus_over_time(graph_HAED, 5, D_RT[:, (3 + (i - 1)):(2 + (i - 1) + 5)], xi_RT, 0.25, gen_com_HA, gen_data, gen_data_conv, offset = (i - 1), HA = true)
        set_optimizer(graph_HAED, Gurobi.Optimizer)
        add_subgraph!(graph, graph_HAED)
    end

    for i in 1:8
        graph_ST = getsubgraphs(graph)[1 + i]
        if i == 1
            ST_sol_time = @elapsed begin
            ST_subgraph_set = getsubgraphs(graph)[2:(i + 1)]
            link_between_ST_sol(graph, graph_ST, ST_subgraph_set, 1; day_num = 1)
            optimize!(graph_ST)
            STUC_sols = value.(all_variables(graph_ST))
            writedlm((@__DIR__)*"/results/32d_serial/STUC_results$(i)_day1.csv", STUC_sols, ',')
            end

            for j in 1:12
                println()
                println("RUNNING $j")
                println()
                link_sol_over_HA(graph, j)
                HA_graph = getsubgraphs(graph)[9 + j]
                optimize!(HA_graph)
                HA_sols = value.(all_variables(HA_graph))
                writedlm((@__DIR__)*"/results/32d_serial/HAED_results$(j)_day1.csv", HA_sols, ',')
            end
        else
            ST_sol_time = @elapsed begin
            ST_subgraph_set = getsubgraphs(graph)[2:(i + 1)]
            link_between_ST_sol(graph, graph_ST, ST_subgraph_set, i; day_num = 1)
            optimize!(graph_ST)
            STUC_sols = value.(all_variables(graph_ST))
            writedlm((@__DIR__)*"/results/32d_serial/STUC_results$(i)_day1.csv", STUC_sols, ',')
            end

            j_range = (1 + (i - 1) * 12):(12 + (i - 1) * 12)
            for j in j_range
                println()
                println("RUNNING $j")
                println()
                link_sol_over_HA(graph, j)
                HA_graph = getsubgraphs(graph)[9 + j]
                optimize!(HA_graph)
                HA_sols = value.(all_variables(HA_graph))
                writedlm((@__DIR__)*"/results/32d_serial/HAED_results$(j)_day1.csv", HA_sols, ',')
            end
        end
    end
    push!(graph_set, graph)
end

function run_day2_serial(graph_set)
    ########################## Solve Day 2 ##################################

    graph_DAUC = OptiGraph()

    build_bus_over_time(graph_DAUC, 25, D_DA[:, 27:51], xi_DA, 1, gen_DA, gen_data_DA, gen_DA, offset = 24)
    link_DA_over_time(graph_DAUC, gen_data_DA)
    graph_DAUC_set = []
    push!(graph_DAUC_set, getsubgraphs(graph_set[1])[1])
    push!(graph_DAUC_set, graph_DAUC)
    println(graph_DAUC_set)
    HA_subgraph_for_DA = getsubgraphs(graph_set[1])[105]
    link_between_DA(graph_DAUC, graph_DAUC_set, 2; HA_graph = HA_subgraph_for_DA)

    set_optimizer(graph_DAUC, Gurobi.Optimizer)

    set_optimizer_attribute(graph_DAUC, "MIPGap", .001)

    println("Optimizing DA Problem")
    optimize!(graph_DAUC)

    graph = OptiGraph()

    add_subgraph!(graph, graph_DAUC)
    set_optimizer(graph, Gurobi.Optimizer)

    DAUC_sols = value.(all_variables(graph_DAUC))
    writedlm((@__DIR__)*"/results/32d_serial/DAUC_results_day2.csv", DAUC_sols, ',')

    for i in 1:8
        local graph_STUC = OptiGraph()
        build_bus_over_time(graph_STUC, 16, D_ST[:, (3 + 96 + (i - 1) * 12):(2 + 96 + (i - 1) * 12 + 16)], xi_ST, 0.25, gen_ST, gen_data, gen_ST, offset = (i - 1) * 12 + 96)
        link_ST_over_time(graph_STUC)
        link_DA_sol_to_ST(graph, graph_STUC, (1 + (i - 1) * 3):(4 + (i - 1) * 3))
        set_optimizer(graph_STUC, Gurobi.Optimizer)
        set_optimizer_attribute(graph_STUC, "MIPGap", .003)
        add_subgraph!(graph, graph_STUC)
    end

    for i in 1:96
        local graph_HAED = OptiGraph()
        build_bus_over_time(graph_HAED, 5, D_RT[:, (96 + 3 + (i - 1)):(96 + 2 + (i - 1) + 5)], xi_RT, 0.25, gen_com_HA, gen_data, gen_data_conv, offset = (i - 1 + 96), HA = true)
        set_optimizer(graph_HAED, Gurobi.Optimizer)
        add_subgraph!(graph, graph_HAED)
    end


    for i in 1:8
        graph_ST = getsubgraphs(graph)[1 + i]

        if i == 1
            subgraph_set1 = getsubgraphs(graph_set[1])[7:9]
            subgraph_set2 = getsubgraphs(graph)[2:(i + 1)]
            ST_subgraph_set = vcat(subgraph_set1, subgraph_set2)
            link_between_ST_sol(graph, graph_ST, ST_subgraph_set, i; day_num = 2, graph_last = graph_set[1])
            optimize!(graph_ST)
            STUC_sols = value.(all_variables(graph_ST))
            writedlm((@__DIR__)*"/results/32d_serial/STUC_results$(i)_day2.csv", STUC_sols, ',')
            link_sol_over_HA_new_day(graph, graph_set[1], monolithic = false)

        elseif i == 2
            subgraph_set1 = getsubgraphs(graph_set[1])[8:9]
            subgraph_set2 = getsubgraphs(graph)[2:(i + 1)]
            ST_subgraph_set = vcat(subgraph_set1, subgraph_set2)
            link_between_ST_sol(graph, graph_ST, ST_subgraph_set, i; day_num = 2, graph_last = graph_set[1])
            optimize!(graph_ST)
            STUC_sols = value.(all_variables(graph_ST))
            writedlm((@__DIR__)*"/results/32d_serial/STUC_results$(i)_day2.csv", STUC_sols, ',')
        elseif i == 3
            subgraph_set1 = getsubgraphs(graph_set[1])[9:9]
            subgraph_set2 = getsubgraphs(graph)[2:(i + 1)]
            ST_subgraph_set = vcat(subgraph_set1, subgraph_set2)
            link_between_ST_sol(graph, graph_ST, ST_subgraph_set, i; day_num = 2, graph_last = graph_set[1])
            optimize!(graph_ST)
            STUC_sols = value.(all_variables(graph_ST))
            writedlm((@__DIR__)*"/results/32d_serial/STUC_results$(i)_day2.csv", STUC_sols, ',')
        else
            ST_subgraph_set = getsubgraphs(graph)[2:(i + 1)]
            link_between_ST_sol(graph, graph_ST, ST_subgraph_set, i; day_num = 2, graph_last = graph_set[1])
            optimize!(graph_ST)
            STUC_sols = value.(all_variables(graph_ST))
            writedlm((@__DIR__)*"/results/32d_serial/STUC_results$(i)_day2.csv", STUC_sols, ',')
        end
        j_range = (1 + (i - 1) * 12):(12 + (i - 1) * 12)
        for j in j_range
            println()
            println("RUNNING $j")
            println()
            link_sol_over_HA(graph, j)
            HA_graph = getsubgraphs(graph)[9 + j]
            optimize!(HA_graph)
            HA_sols = value.(all_variables(HA_graph))
            writedlm((@__DIR__)*"/results/32d_serial/HAED_results$(j)_day2.csv", HA_sols, ',')
        end
    end

    push!(graph_set, graph)
end

function run_dayi_serial(graph_set, day)
    if length(graph_set) >= 3
        popfirst!(graph_set)
    end
    ########################## Solve Day i ##################################

    graph_DAUC = OptiGraph()

    DA_range = (3 + 24 * (day - 1)):(3 + 24 * day)
    build_bus_over_time(graph_DAUC, 25, D_DA[:, DA_range], xi_DA, 1, gen_DA, gen_data_DA, gen_DA, offset = (day - 1) * 24)
    link_DA_over_time(graph_DAUC, gen_data_DA)
    graph_DAUC_set = []
    for k in 1:length(graph_set)
        push!(graph_DAUC_set, getsubgraphs(graph_set[k])[1])
    end
    push!(graph_DAUC_set, graph_DAUC)
    println(graph_DAUC_set)
    HA_subgraph_for_DA = getsubgraphs(graph_set[length(graph_set)])[105]
    link_between_DA(graph_DAUC, graph_DAUC_set, day; HA_graph = HA_subgraph_for_DA)

    set_optimizer(graph_DAUC, Gurobi.Optimizer)

    if day == 27 || day == 3 || day == 31
        set_optimizer_attribute(graph_DAUC, "MIPGap", .0005)
    elseif day == 26
        set_optimizer_attribute(graph_DAUC, "MIPGap", .003)
    else
        set_optimizer_attribute(graph_DAUC, "MIPGap", .001)
    end
    set_optimizer_attribute(graph_DAUC, "Threads", 8)

    println("Optimizing DA Problem")
    optimize!(graph_DAUC)

    graph = OptiGraph()

    add_subgraph!(graph, graph_DAUC)
    set_optimizer(graph, Gurobi.Optimizer)

    DAUC_sols = value.(all_variables(graph_DAUC))
    writedlm((@__DIR__)*"/results/32d_serial/DAUC_results_day$(day).csv", DAUC_sols, ',')

    for i in 1:8
        local graph_STUC = OptiGraph()
        ST_range = (3 + (day - 1) * 96 + (i - 1) * 12):(2 + (day - 1) * 96 + (i - 1) * 12 + 16)
        build_bus_over_time(graph_STUC, 16, D_ST[:, ST_range], xi_ST, 0.25, gen_ST, gen_data, gen_ST, offset = (day - 1) * 96 + (i - 1) * 12)
        link_ST_over_time(graph_STUC)
        link_DA_sol_to_ST(graph, graph_STUC, (1 + (i - 1) * 3):(4 + (i - 1) * 3))
        set_optimizer(graph_STUC, Gurobi.Optimizer)
        if (day == 28) || ((day == 31) && (i == 2))
            set_optimizer_attribute(graph_STUC, "MIPGap", .001)
        elseif day == 26#day == 27 ||
            set_optimizer_attribute(graph_STUC, "MIPGap", .0005)
        elseif ((day == 27) && (i == 6)) || ((day == 31) && (i == 1))
            set_optimizer_attribute(graph_STUC, "MIPGap", .004)
        elseif (day == 6) && (i == 3)
            set_optimizer_attribute(graph_STUC, "MIPGap", .006)
        else
            set_optimizer_attribute(graph_STUC, "MIPGap", .003)
        end
        set_optimizer_attribute(graph_STUC, "Threads", 8)

        add_subgraph!(graph, graph_STUC)
    end

    for i in 1:96
        local graph_HAED = OptiGraph()
        HA_range = ((day - 1) * 96 + 3 + (i - 1)):((day - 1) * 96 + 2 + (i - 1) + 5)
        build_bus_over_time(graph_HAED, 5, D_RT[:, HA_range], xi_RT, 0.25, gen_com_HA, gen_data, gen_data_conv, offset = ((day - 1) * 96 + (i - 1)), HA = true)
        set_optimizer(graph_HAED, Gurobi.Optimizer)
        add_subgraph!(graph, graph_HAED)
    end


    for i in 1:8
        graph_ST = getsubgraphs(graph)[1 + i]

        if i == 1
            subgraph_set1 = getsubgraphs(graph_set[2])[7:9]
            subgraph_set2 = getsubgraphs(graph)[2:(i + 1)]
            ST_subgraph_set = vcat(subgraph_set1, subgraph_set2)
            link_between_ST_sol(graph, graph_ST, ST_subgraph_set, i; day_num = day, graph_last = graph_set[2])
            optimize!(graph_ST)
            STUC_sols = value.(all_variables(graph_ST))
            writedlm((@__DIR__)*"/results/32d_serial/STUC_results$(i)_day$(day).csv", STUC_sols, ',')
            link_sol_over_HA_new_day(graph, graph_set[2], monolithic = false)

        elseif i == 2
            subgraph_set1 = getsubgraphs(graph_set[2])[8:9]
            subgraph_set2 = getsubgraphs(graph)[2:(i + 1)]
            ST_subgraph_set = vcat(subgraph_set1, subgraph_set2)
            link_between_ST_sol(graph, graph_ST, ST_subgraph_set, i; day_num = day, graph_last = graph_set[2])
            optimize!(graph_ST)
            STUC_sols = value.(all_variables(graph_ST))
            writedlm((@__DIR__)*"/results/32d_serial/STUC_results$(i)_day$(day).csv", STUC_sols, ',')
        elseif i == 3
            subgraph_set1 = getsubgraphs(graph_set[2])[9:9]
            subgraph_set2 = getsubgraphs(graph)[2:(i + 1)]
            ST_subgraph_set = vcat(subgraph_set1, subgraph_set2)
            link_between_ST_sol(graph, graph_ST, ST_subgraph_set, i; day_num = day, graph_last = graph_set[2])
            optimize!(graph_ST)
            STUC_sols = value.(all_variables(graph_ST))
            writedlm((@__DIR__)*"/results/32d_serial/STUC_results$(i)_day$(day).csv", STUC_sols, ',')
        else
            ST_subgraph_set = getsubgraphs(graph)[2:(i + 1)]
            link_between_ST_sol(graph, graph_ST, ST_subgraph_set, i; day_num = day, graph_last = graph_set[2])
            optimize!(graph_ST)
            STUC_sols = value.(all_variables(graph_ST))
            writedlm((@__DIR__)*"/results/32d_serial/STUC_results$(i)_day$(day).csv", STUC_sols, ',')
        end
        j_range = (1 + (i - 1) * 12):(12 + (i - 1) * 12)
        for j in j_range
            println()
            println("RUNNING $j")
            println()
            link_sol_over_HA(graph, j)
            HA_graph = getsubgraphs(graph)[9 + j]
            optimize!(HA_graph)
            HA_sols = value.(all_variables(HA_graph))
            writedlm((@__DIR__)*"/results/32d_serial/HAED_results$(j)_day$(day).csv", HA_sols, ',')
        end
    end

    push!(graph_set, graph)
end

graph_set = []
d1_time = @elapsed run_day1_serial(graph_set)
d2_time = @elapsed run_day2_serial(graph_set)

for i in 3:32
    di_time = @elapsed run_dayi_serial(graph_set, i)
end
