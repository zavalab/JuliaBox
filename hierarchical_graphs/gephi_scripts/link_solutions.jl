using Plasmo
using DelimitedFiles, CSV, DataFrames
using Gurobi


function link_DA_sol_to_ST(graph::OptiGraph, ST_subgraph, tau_0 = 1:4, dt = .25, gen_com_set = gen_ST)

    DA_subgraph = getsubgraphs(graph)[1]
    #ST_subgraph = getsubgraphs(graph)[1 + ST_sub_num]

    num_DA_gen = size(gen_DA, 1)

    for (i, index) in enumerate(tau_0)
        if index <= length(getsubgraphs(DA_subgraph)) # If index is less than the number of time points in DAUC problem
            DA_subgraph_i = getsubgraphs(DA_subgraph)[index]
            ST_subgraphs_i = getsubgraphs(ST_subgraph)[(1 + (i - 1) * 4):(i * 4)]

            for k in 1:num_DA_gen
                gen_name = gen_DA[k, "Generator Name"]
                bus_name = gen_DA[k, "Bus of Connection"]
                max_ramp_up = gen_DA[k, "Max Ramp Up (MW/min)"] * 60 * dt
                max_ramp_down = gen_DA[k, "Max Ramp Down (MW/min)"] * 60 * dt

                max_startup = gen_DA[k, "Max Capacity (MW)"]
                max_shutdown = gen_DA[k, "Max Capacity (MW)"]

                gen_max_cap = gen_DA[k, "Max Capacity (MW)"]
                gen_min_cap = gen_DA[k, "Min Stable Level (MW)"]

                for l in 1:4
                    @linkconstraint(graph, (ST_subgraphs_i[l][:bus][bus_name][:Gen_p][gen_name] + ST_subgraphs_i[l][:bus][bus_name][:Gen_m][gen_name])
                                            - (value(DA_subgraph_i[:bus][bus_name][:Gen_p][gen_name]) + value(DA_subgraph_i[:bus][bus_name][:Gen_m][gen_name])) <= max_ramp_up / dt
                    )
                    @linkconstraint(graph, - (ST_subgraphs_i[l][:bus][bus_name][:Gen_p][gen_name] + ST_subgraphs_i[l][:bus][bus_name][:Gen_m][gen_name])
                                            + (value(DA_subgraph_i[:bus][bus_name][:Gen_p][gen_name]) + value(DA_subgraph_i[:bus][bus_name][:Gen_m][gen_name])) <= max_ramp_up / dt
                    )
                end


                for l in 2:4
                    bus_next = ST_subgraphs_i[l][:bus][bus_name]
                    bus_last = ST_subgraphs_i[l - 1][:bus][bus_name]
                    @linkconstraint(ST_subgraph, (bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name]) -
                                           (bus_last[:Gen_p][gen_name] + bus_last[:Gen_m][gen_name]) <=
                                           max_ramp_up
                    )
                    @linkconstraint(ST_subgraph, (bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name]) -
                                                 (bus_last[:Gen_p][gen_name] + bus_last[:Gen_m][gen_name]) >=
                                                  - max_ramp_down
                    )
                end

                if i != 1
                    bus_next = getsubgraphs(ST_subgraph)[1 + (i - 1) * 4][:bus][bus_name]
                    bus_last = getsubgraphs(ST_subgraph)[(i - 1) * 4][:bus][bus_name]
                    bus_next_DA = getsubgraphs(DA_subgraph)[index][:bus][bus_name]
                    bus_last_DA = getsubgraphs(DA_subgraph)[index - 1][:bus][bus_name]

                    @linkconstraint(ST_subgraph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                        bus_last[:Gen_p][gen_name] - bus_last[:Gen_m][gen_name] <=
                        (max_startup - max_ramp_up - gen_min_cap) * value(bus_next_DA[:s][gen_name]) +
                        (max_ramp_up + gen_min_cap) * value(bus_next_DA[:x][gen_name]) - gen_min_cap * value(bus_last_DA[:x][gen_name])
                    )

                    @linkconstraint(ST_subgraph, bus_last[:Gen_p][gen_name] + bus_last[:Gen_m][gen_name] -
                        bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                        (max_shutdown - max_ramp_down - gen_min_cap) * value(bus_next_DA[:z][gen_name]) +
                        (max_ramp_down + gen_min_cap) * value(bus_last_DA[:x][gen_name]) - gen_min_cap * value(bus_next_DA[:x][gen_name])
                    )
                end
                # TO CHECK: max_ramp down or up should maybe be divided by dt;
                for (j, subgraph) in enumerate(ST_subgraphs_i)
                    @linkconstraint(ST_subgraph, gen_max_cap * value(DA_subgraph_i[:bus][bus_name][:x][gen_name]) >=
                                           subgraph[:bus][bus_name][:Gen_p][gen_name] + subgraph[:bus][bus_name][:Gen_m][gen_name]
                    )
                    @linkconstraint(ST_subgraph, gen_min_cap * value(DA_subgraph_i[:bus][bus_name][:x][gen_name]) <=
                                           subgraph[:bus][bus_name][:Gen_p][gen_name] + subgraph[:bus][bus_name][:Gen_m][gen_name]
                    )
                end
            end
        end
    end
end

function link_between_ST_sol(graph, ST_subgraph, ST_graphs, i; dt = .25, day_num = 1, graph_last = nothing)
    # ST_subgraph is the STUC_graph for a 4 hour time interval
    # ST_graphs is the set of STUC_graphs for time 1 to time i

    function get_s_list(ST_graphs, bus_name, gen_name, i)
        if i == 1
            s_list = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1]
        elseif i == 2
            s_list1 = Vector{Any}([getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1])
            s_list2 = Vector{Any}([value(getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:s][gen_name]) for j in 12:-1:1])
            s_list = vcat(s_list1, s_list2)
        elseif i == 3
            s_list1 = Vector{Any}([getsubgraphs(ST_graphs[3])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1])
            s_list2 = Vector{Any}([value(getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:s][gen_name]) for j in 12:-1:1])
            s_list3 = Vector{Any}([value(getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:s][gen_name]) for j in 12:-1:1])
            s_list = vcat(s_list1, s_list2, s_list3)
        else
            s_list1 = Vector{Any}([getsubgraphs(ST_graphs[i])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1])
            s_list2 = Vector{Any}([value(getsubgraphs(ST_graphs[i - 1])[j][:bus][bus_name][:s][gen_name]) for j in 12:-1:1])
            s_list3 = Vector{Any}([value(getsubgraphs(ST_graphs[i - 2])[j][:bus][bus_name][:s][gen_name]) for j in 12:-1:1])
            s_list4 = Vector{Any}([value(getsubgraphs(ST_graphs[i - 3])[j][:bus][bus_name][:s][gen_name]) for j in 12:-1:1])
            s_list = vcat(s_list1, s_list2, s_list3, s_list4)
        end
        return s_list
    end

    function get_z_list(ST_graphs, bus_name, gen_name, i)
        if i == 1
            z_list = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1]
        elseif i == 2
            z_list1 = Vector{Any}([getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1])
            z_list2 = Vector{Any}([value(getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:z][gen_name]) for j in 12:-1:1])
            z_list = vcat(z_list1, z_list2)
        elseif i == 3
            z_list1 = Vector{Any}([getsubgraphs(ST_graphs[3])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1])
            z_list2 = Vector{Any}([value(getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:z][gen_name]) for j in 12:-1:1])
            z_list3 = Vector{Any}([value(getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:z][gen_name]) for j in 12:-1:1])
            z_list = vcat(z_list1, z_list2, z_list3)
        else
            z_list1 = Vector{Any}([getsubgraphs(ST_graphs[i])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1])
            z_list2 = Vector{Any}([value(getsubgraphs(ST_graphs[i - 1])[j][:bus][bus_name][:z][gen_name]) for j in 12:-1:1])
            z_list3 = Vector{Any}([value(getsubgraphs(ST_graphs[i - 2])[j][:bus][bus_name][:z][gen_name]) for j in 12:-1:1])
            z_list4 = Vector{Any}([value(getsubgraphs(ST_graphs[i - 3])[j][:bus][bus_name][:z][gen_name]) for j in 12:-1:1])
            z_list = vcat(z_list1, z_list2, z_list3, z_list4)
        end
        return z_list
    end

    function link_subgraphs_uptime(graph, ST_graphs, min_up_time, bus_name, gen_name, s_list, i)
        if i == 1
            if min_up_time >= 4
                for j in 1:15
                    subgraph_index = 17 - j
                    s_range = (j):16
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                #@linkconstraint(graph, sum(s_list[2:16]) <= getsubgraphs(ST_graphs)[16][:bus][bus_name][:x][gen_name])
            else
                up_time_diff = Int(4 * (4 - min_up_time))
                for j in 1:up_time_diff
                    s_range = (j):(Int(min_up_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):15 # this maybe should go to 15, not 16
                    s_range = (j):16
                    subgraph_index = 17 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 2 # s_list is length 28
            if (min_up_time >= 3) && (min_up_time < 7)
                up_time_diff = Int(4 * (7 - min_up_time))

                for j in 1:up_time_diff
                    s_range = (j):(Int(min_up_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):16
                    s_range = (j):28
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            elseif (min_up_time >= 7)
                for j in 1:16
                    subgraph_index = 17 - j
                    s_range = (j):28
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:16
                    subgraph_index = 17 - j
                    s_range = (j):(Int(min_up_time * 4 + j - 1))
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 3 #s_list is length 40
            if (min_up_time >= 6)
                up_time_diff = Int(4 * (10 - min_up_time))
                for j in 1:up_time_diff
                    s_range = (j):(Int(min_up_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):16
                    s_range = (j):40
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:16
                    subgraph_index = 17 - j
                    s_range = (j):(Int(min_up_time * 4 + j - 1))
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        else #s_list is length 52
            for j in 1:16
                subgraph_index = 17 - j
                s_range = (j):(Int(min_up_time * 4 + j - 1))
                @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[i])[subgraph_index][:bus][bus_name][:x][gen_name])
            end
        end
    end

    function link_subgraphs_downtime(graph, ST_graphs, min_down_time, bus_name, gen_name, z_list, i)
        if i == 1
            if min_down_time >= 4
                for j in 1:15
                    subgraph_index = 17 - j
                    z_range = j:16
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                down_time_diff = Int(4 * (4 - min_down_time))
                for j in 1:down_time_diff
                    z_range = j:(Int(min_down_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (down_time_diff + 1):16
                    z_range = j:16
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 2 # s_list is length 28
            if (min_down_time >= 3) && (min_down_time < 7)
                down_time_diff = Int(4 * (7 - min_down_time))

                for j in 1:down_time_diff
                    z_range = j:(Int(min_down_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (down_time_diff + 1):16
                    z_range = j:28
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            elseif (min_down_time >= 7)
                for j in 1:16
                    subgraph_index = 17 - j
                    z_range = j:28
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                #@linkconstraint(graph, sum(z_list[1:28]) <= getsubgraphs(ST_graphs[2])[16][:bus][bus_name][:x][gen_name])
            else
                for j in 1:16
                    subgraph_index = 17 - j
                    z_range = j:(Int(min_down_time * 4 + j - 1))
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 3 #s_list is length 40
            if (min_down_time >= 6)
                down_time_diff = Int(4 * (10 - min_down_time))
                for j in 1:down_time_diff
                    z_range = j:(Int(min_down_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (down_time_diff + 1):16
                    z_range = j:40
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:16
                    subgraph_index = 17 - j
                    z_range = j:(Int(min_down_time * 4 + j - 1))
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        else #s_list is length 52
            for j in 1:16
                subgraph_index = 17 - j
                z_range = j:(Int(min_down_time * 4 + j - 1))
                @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(ST_graphs[i])[subgraph_index][:bus][bus_name][:x][gen_name])
            end
        end
    end

    ST_graphs_len = length(ST_graphs)
    num_ST_gen = size(gen_ST, 1)

    # Need a constraint added for the DA generators; if i != 1, then there is a constraint between ST subproblem 1 and 2, and 2 and 3, and so forth;
    for j in 1:num_ST_gen
        gen_name = gen_ST[j, "Generator Name"]
        bus_name = gen_ST[j, "Bus of Connection"]

        min_up_time = gen_ST[j, "Min Up Time (h)"]
        min_down_time = gen_ST[j, "Min Down Time (h)"]

        max_ramp_up = gen_ST[j, "Max Ramp Up (MW/min)"] * 60 * dt
        max_ramp_down = gen_ST[j, "Max Ramp Down (MW/min)"] * 60 * dt

        max_cap = gen_ST[j, "Max Capacity (MW)"]
        min_cap = gen_ST[j, "Min Stable Level (MW)"]

        max_startup = gen_ST[j, "Max Capacity (MW)"]
        max_shutdown = gen_ST[j, "Max Capacity (MW)"]

        s_list = get_s_list(ST_graphs, bus_name, gen_name, ST_graphs_len)
        z_list = get_z_list(ST_graphs, bus_name, gen_name, ST_graphs_len)
        if !(isequal(0, min_up_time))
            link_subgraphs_uptime(ST_subgraph, ST_graphs, min_up_time, bus_name, gen_name, s_list, ST_graphs_len)
        end
        if !(isequal(0, min_down_time))
            link_subgraphs_downtime(ST_subgraph, ST_graphs, min_down_time, bus_name, gen_name, z_list, ST_graphs_len)
        end

        if ((day_num != 1) || (i != 1))

            if i == 1
                bus_next = getsubgraphs(ST_subgraph)[1][:bus][bus_name]
                bus_last = getsubgraphs(getsubgraphs(graph_last)[9])[12][:bus][bus_name]
                HA_subgraph = getsubgraphs(graph_last)[105]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]
            else
                bus_next = getsubgraphs(ST_graphs[ST_graphs_len])[1][:bus][bus_name]
                bus_last = getsubgraphs(ST_graphs[ST_graphs_len - 1])[12][:bus][bus_name]
                HA_subgraph = getsubgraphs(graph)[9 + (i - 1) * 12]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]
            end
            @linkconstraint(ST_subgraph, bus_next[:x][gen_name] - value(bus_last[:x][gen_name]) == bus_next[:s][gen_name] - bus_next[:z][gen_name])

            @linkconstraint(ST_subgraph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                value(bus_last_HA[:Gen_p][gen_name]) - value(bus_last_HA[:Gen_m][gen_name]) <=
                (max_startup - max_ramp_up - min_cap) * bus_next[:s][gen_name] +
                (max_ramp_up + min_cap) * bus_next[:x][gen_name] - min_cap * value(bus_last[:x][gen_name])
            )

            @linkconstraint(ST_subgraph, value(bus_last_HA[:Gen_p][gen_name]) + value(bus_last_HA[:Gen_m][gen_name]) -
                bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                (max_shutdown - max_ramp_down - min_cap) * bus_next[:z][gen_name] +
                (max_ramp_down + min_cap) * value(bus_last[:x][gen_name]) - min_cap * bus_next[:x][gen_name]
            )
        end
    end

    num_DA_gen = size(gen_DA, 1)

    if (i != 1) || (day_num != 1)
        for j in 1:num_DA_gen
            gen_name = gen_DA[j, "Generator Name"]
            bus_name = gen_DA[j, "Bus of Connection"]

            min_up_time = gen_DA[j, "Min Up Time (h)"]
            min_down_time = gen_DA[j, "Min Down Time (h)"]

            max_ramp_up = gen_DA[j, "Max Ramp Up (MW/min)"] * 60 * dt
            max_ramp_down = gen_DA[j, "Max Ramp Down (MW/min)"] * 60 * dt

            max_cap = gen_DA[j, "Max Capacity (MW)"]
            min_cap = gen_DA[j, "Min Stable Level (MW)"]

            max_startup = gen_DA[j, "Max Capacity (MW)"]
            max_shutdown = gen_DA[j, "Max Capacity (MW)"]

            if i == 1
                bus_next = getsubgraphs(ST_subgraph)[1][:bus][bus_name]
                bus_last = getsubgraphs(getsubgraphs(graph_last)[9])[12][:bus][bus_name]

                HA_subgraph = getsubgraphs(graph_last)[105]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]

                DA_subgraph_last = getsubgraphs(graph_last)[1]
                bus_last_DA = getsubgraphs(DA_subgraph_last)[24][:bus][bus_name]

                DA_subgraph_next = getsubgraphs(graph)[1]
                bus_next_DA = getsubgraphs(DA_subgraph_next)[1][:bus][bus_name]
            else
                bus_next = getsubgraphs(ST_subgraph)[1][:bus][bus_name]
                bus_last = getsubgraphs(ST_graphs[ST_graphs_len - 1])[12][:bus][bus_name]

                HA_subgraph = getsubgraphs(graph)[9 + (i - 1) * 12]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]

                DA_subgraph = getsubgraphs(graph)[1]
                bus_last_DA = getsubgraphs(DA_subgraph)[(i - 1) * 3][:bus][bus_name]

                DA_subgraph = getsubgraphs(graph)[1]
                bus_next_DA = getsubgraphs(DA_subgraph)[1 + (i - 1) * 3][:bus][bus_name]
            end

            @linkconstraint(ST_subgraph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                value(bus_last_HA[:Gen_p][gen_name]) - value(bus_last_HA[:Gen_m][gen_name]) <=
                (max_startup - max_ramp_up / dt - min_cap) * value(bus_next_DA[:s][gen_name]) +
                (max_ramp_up / dt + min_cap) * value(bus_next_DA[:x][gen_name]) - min_cap * value(bus_last_DA[:x][gen_name])
            )

            @linkconstraint(ST_subgraph, value(bus_last_HA[:Gen_p][gen_name]) + value(bus_last_HA[:Gen_m][gen_name]) -
                bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                (max_shutdown - max_ramp_down / dt - min_cap) * value(bus_next_DA[:z][gen_name]) +
                (max_ramp_down / dt + min_cap) * value(bus_last_DA[:x][gen_name]) - min_cap * value(bus_next_DA[:x][gen_name])
            )
        end
    end
end

function link_sol_over_HA(graph, tau = 1; dt = .25)

    ST_layer = Int(ceil(tau / 12))
    ST_start_index = tau%12

    if ST_start_index == 0
        ST_start_index = 12
    end

    DA_subgraph = getsubgraphs(graph)[1]
    ST_subgraph = getsubgraphs(graph)[1 + ST_layer]
    HA_subgraph = getsubgraphs(graph)[9 + tau]
    tau_range = ST_start_index:(ST_start_index + 4)
    for (i, index) in enumerate(tau_range)
        DA_index = Int(ceil((tau + i - 1) / 4))

        DA_subgraph_i = getsubgraphs(DA_subgraph)[DA_index]
        ST_subgraph_i = getsubgraphs(ST_subgraph)[index]
        HA_subgraph_i = getsubgraphs(HA_subgraph)[i]

        num_conv_gen = size(gen_data_conv, 1)

        DA_len = size(gen_DA, 1)
        for k in 1:DA_len
            gen_name = gen_DA[k, "Generator Name"]
            bus_name = gen_DA[k, "Bus of Connection"]

            gen_max_cap = gen_DA[k, "Max Capacity (MW)"]
            gen_min_cap = gen_DA[k, "Min Stable Level (MW)"]

            max_ramp_up = gen_DA[k, "Max Ramp Up (MW/min)"] * 60 * dt# * dt #TODO: check if this should be multiplied by dt or not
            max_ramp_down = gen_DA[k, "Max Ramp Down (MW/min)"] * 60 * dt# * dt

            max_startup = gen_DA[k, "Max Capacity (MW)"]
            max_shutdown = gen_DA[k, "Max Capacity (MW)"]

            if !(isequal(0, max_ramp_up))
                @linkconstraint(HA_subgraph, (value(DA_subgraph_i[:bus][bus_name][:Gen_p][gen_name]) + value(DA_subgraph_i[:bus][bus_name][:Gen_m][gen_name])) -
                                       (HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up / dt
                )
                @linkconstraint(HA_subgraph, -(value(DA_subgraph_i[:bus][bus_name][:Gen_p][gen_name]) + value(DA_subgraph_i[:bus][bus_name][:Gen_m][gen_name])) +
                                       (HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up / dt
                )
            end

            @linkconstraint(HA_subgraph, gen_max_cap * value(DA_subgraph_i[:bus][bus_name][:x][gen_name]) >=
                                   HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]
            )
            @linkconstraint(HA_subgraph, gen_min_cap * value(DA_subgraph_i[:bus][bus_name][:x][gen_name]) <=
                                   HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]
            )

            if (tau != 1) && (DA_index != 1)
                bus_next = HA_subgraph_i[:bus][bus_name]
                if i != 1
                    bus_last = getsubgraphs(HA_subgraph)[i - 1][:bus][bus_name]
                    bus_last_gen_p = bus_last[:Gen_p][gen_name]
                    bus_last_gen_m = bus_last[:Gen_m][gen_name]
                else
                    HA_subgraph_last = getsubgraphs(graph)[9 + tau - 1]
                    bus_last = getsubgraphs(HA_subgraph_last)[1][:bus][bus_name]
                    bus_last_gen_p = value(bus_last[:Gen_p][gen_name])
                    bus_last_gen_m = value(bus_last[:Gen_m][gen_name])
                end

                if (tau + i - 1) % 4 == 1 ### THIS IS THE SECTION THAT NEEDS TO BE CHANGED TO AVOID INFEASIBILITY
                    bus_next_DA = getsubgraphs(DA_subgraph)[DA_index][:bus][bus_name]
                    bus_last_DA = getsubgraphs(DA_subgraph)[DA_index - 1][:bus][bus_name]
                    @linkconstraint(HA_subgraph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                        bus_last_gen_p - bus_last_gen_m <=
                        (max_startup - max_ramp_up - gen_min_cap) * value(bus_next_DA[:s][gen_name]) +
                        (max_ramp_up + gen_min_cap) * value(bus_next_DA[:x][gen_name]) - gen_min_cap * value(bus_last_DA[:x][gen_name])
                    )

                    @linkconstraint(HA_subgraph, bus_last_gen_p + bus_last_gen_m -
                        bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                        (max_shutdown - max_ramp_down - gen_min_cap) * value(bus_next_DA[:z][gen_name]) +
                        (max_ramp_down + gen_min_cap) * value(bus_last_DA[:x][gen_name]) - gen_min_cap * value(bus_next_DA[:x][gen_name])
                    )
                else
                    @linkconstraint(HA_subgraph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                        bus_last_gen_p - bus_last_gen_m <= max_ramp_up
                    )
                    @linkconstraint(HA_subgraph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                        bus_last_gen_p - bus_last_gen_m >= - max_ramp_down
                    )
                end
            end
        end

        ST_len = size(gen_ST, 1)
        for k in 1:ST_len
            gen_name = gen_ST[k, "Generator Name"]
            bus_name = gen_ST[k, "Bus of Connection"]

            gen_max_cap = gen_ST[k, "Max Capacity (MW)"]
            gen_min_cap = gen_ST[k, "Min Stable Level (MW)"]

            max_ramp_up = gen_ST[k, "Max Ramp Up (MW/min)"] * 60 * dt
            max_ramp_down = gen_ST[k, "Max Ramp Down (MW/min)"] * 60 * dt

            max_startup = gen_ST[k, "Max Capacity (MW)"]
            max_shutdown = gen_ST[k, "Max Capacity (MW)"]

            if !(isequal(0, max_ramp_up))
                @linkconstraint(HA_subgraph, (value(ST_subgraph_i[:bus][bus_name][:Gen_p][gen_name]) + value(ST_subgraph_i[:bus][bus_name][:Gen_m][gen_name])) -
                                       (HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up# / dt
                )
                @linkconstraint(HA_subgraph, -(value(ST_subgraph_i[:bus][bus_name][:Gen_p][gen_name]) + value(ST_subgraph_i[:bus][bus_name][:Gen_m][gen_name])) +
                                       (HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up# / dt
                )
            end

            @linkconstraint(HA_subgraph, gen_max_cap * value(ST_subgraph_i[:bus][bus_name][:x][gen_name]) >=
                                   HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]
            )
            @linkconstraint(HA_subgraph, gen_min_cap * value(ST_subgraph_i[:bus][bus_name][:x][gen_name]) <=
                                   HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]
            )

            if (tau != 1)
                bus_next = HA_subgraph_i[:bus][bus_name]
                if i != 1
                    bus_last = getsubgraphs(HA_subgraph)[i - 1][:bus][bus_name]
                    bus_last_gen_p = bus_last[:Gen_p][gen_name]
                    bus_last_gen_m = bus_last[:Gen_m][gen_name]
                else
                    HA_subgraph_last = getsubgraphs(graph)[9 + tau - 1]
                    bus_last = getsubgraphs(HA_subgraph_last)[1][:bus][bus_name]
                    bus_last_gen_p = value(bus_last[:Gen_p][gen_name])
                    bus_last_gen_m = value(bus_last[:Gen_m][gen_name])
                end
                bus_next_ST = getsubgraphs(ST_subgraph)[index][:bus][bus_name]
                if index != 1 # Previously was ST_start_index != 1
                    bus_last_ST = getsubgraphs(ST_subgraph)[index - 1][:bus][bus_name]
                else
                    ST_subgraph_last = getsubgraphs(graph)[1 + ST_layer - 1]
                    bus_last_ST = getsubgraphs(ST_subgraph_last)[12][:bus][bus_name]
                end
                @linkconstraint(HA_subgraph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                    bus_last_gen_p - bus_last_gen_m <=
                    (max_startup - max_ramp_up - gen_min_cap) * value(bus_next_ST[:s][gen_name]) +
                    (max_ramp_up + gen_min_cap) * value(bus_next_ST[:x][gen_name]) - gen_min_cap * value(bus_last_ST[:x][gen_name])
                )

                @linkconstraint(HA_subgraph, bus_last_gen_p + bus_last_gen_m -
                    bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                    (max_shutdown - max_ramp_down - gen_min_cap) * value(bus_next_ST[:z][gen_name]) +
                    (max_ramp_down + gen_min_cap) * value(bus_last_ST[:x][gen_name]) - gen_min_cap * value(bus_next_ST[:x][gen_name])
                )
            elseif i != 1
                bus_next = HA_subgraph_i[:bus][bus_name]
                bus_last = getsubgraphs(HA_subgraph)[i - 1][:bus][bus_name]
                bus_last_gen_p = bus_last[:Gen_p][gen_name]
                bus_last_gen_m = bus_last[:Gen_m][gen_name]

                bus_next_ST = getsubgraphs(ST_subgraph)[index][:bus][bus_name]
                bus_last_ST = getsubgraphs(ST_subgraph)[index - 1][:bus][bus_name]
                @linkconstraint(HA_subgraph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                    bus_last_gen_p - bus_last_gen_m <=
                    (max_startup - max_ramp_up - gen_min_cap) * value(bus_next_ST[:s][gen_name]) +
                    (max_ramp_up + gen_min_cap) * value(bus_next_ST[:x][gen_name]) - gen_min_cap * value(bus_last_ST[:x][gen_name])
                )

                @linkconstraint(HA_subgraph, bus_last_gen_p + bus_last_gen_m -
                    bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                    (max_shutdown - max_ramp_down - gen_min_cap) * value(bus_next_ST[:z][gen_name]) +
                    (max_ramp_down + gen_min_cap) * value(bus_last_ST[:x][gen_name]) - gen_min_cap * value(bus_next_ST[:x][gen_name])
                )
            end
        end
    end
end


function link_between_DA(DA_graph, DA_graph_set, day_num; dt = 1, HA_graph = nothing) # graph here is the DAUC_graph; DA_graph_set is a length 2 or length 3 list of DAUC_graphs

    function get_s_list(DA_graph_set, bus_name, gen_name, day_num)

        if day_num == 1
            s_list = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:s][gen_name] for j in 25:-1:1]
        elseif day_num == 2
            s_list1 = [getsubgraphs(DA_graph_set[2])[j][:bus][bus_name][:s][gen_name] for j in 25:-1:1]
            s_list2 = [value(getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:s][gen_name]) for j in 24:-1:1]
            s_list = vcat(s_list1, s_list2)
        else
            s_list1 = [getsubgraphs(DA_graph_set[3])[j][:bus][bus_name][:s][gen_name] for j in 25:-1:1]
            s_list2 = [value(getsubgraphs(DA_graph_set[2])[j][:bus][bus_name][:s][gen_name]) for j in 24:-1:1]
            s_list3 = [value(getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:s][gen_name]) for j in 24:-1:1]
            s_list = vcat(s_list1, s_list2, s_list3)
        end
        return s_list
    end

    function get_z_list(DA_graph_set, bus_name, gen_name, day_num)
        if day_num == 1
            z_list = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:z][gen_name] for j in 25:-1:1]
        elseif day_num == 2
            z_list1 = [getsubgraphs(DA_graph_set[2])[j][:bus][bus_name][:z][gen_name] for j in 25:-1:1]
            z_list2 = [value(getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:z][gen_name]) for j in 24:-1:1]
            z_list = vcat(z_list1, z_list2)
        else
            z_list1 = [getsubgraphs(DA_graph_set[3])[j][:bus][bus_name][:z][gen_name] for j in 25:-1:1]
            z_list2 = [value(getsubgraphs(DA_graph_set[2])[j][:bus][bus_name][:z][gen_name]) for j in 24:-1:1]
            z_list3 = [value(getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:z][gen_name]) for j in 24:-1:1]
            z_list = vcat(z_list1, z_list2, z_list3)
        end
        return z_list
    end

    function link_subgraphs_uptime(graph, DA_graphs, min_up_time, bus_name, gen_name, s_list, i)
        if i == 1# Look at code for elseif i == 3 below
            # if min_up_time is > 24, then sum over all available time points up to current time point
            # (e.g., if 36, than time points 24-36 (and 37) will sum over all previous time points, then time point 38 will sum over time points
            # 2:37 and 39 sums over time points 3:38
            # Correction to what is above: you only need to sum over 2:36 for time point 37, then sum over 3:37 for time point 38
            # TODO: I need to double check all my up time summations; I might be summing over one point too many;
            # For example, if I have an up time of 4 hours, at hour 6 I need to sum over hours 3, 4, and 5; I do not need to include 2 because
            # if it turned on in 2 (s = 1), then I could shut off in time point 6
            if min_up_time >24
                for j in 1:24
                    subgraph_index = 26 - j
                    s_range = (j):25
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                #@linkconstraint(graph, sum(s_list[2:16]) <= getsubgraphs(ST_graphs)[16][:bus][bus_name][:x][gen_name])
            else
                up_time_diff = Int(25 - min_up_time)
                for j in 1:up_time_diff
                    s_range = (j):(Int(min_up_time + j - 1))
                    subgraph_index = 26 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):24
                    s_range = (j):25
                    subgraph_index = 26 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 2 # s_list is length 28
            if min_up_time > 24
                up_time_diff = Int(49 - min_up_time)

                for j in 1:up_time_diff
                    s_range = (j):(Int(min_up_time + j - 1))
                    subgraph_index = 26 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):25
                    s_range = (j):49
                    subgraph_index = 26 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:25
                    subgraph_index = 26 - j
                    s_range = (j):(Int(min_up_time + j - 1))
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        else
            for j in 1:25
                subgraph_index = 26 - j
                s_range = (j):(Int(min_up_time + j - 1))
                @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
            end
        end
    end

    function link_subgraphs_downtime(graph, DA_graphs, min_down_time, bus_name, gen_name, z_list, i)
        if i == 1
            if min_down_time > 24
                for j in 1:24
                    subgraph_index = 26 - j
                    z_range = j:25
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                down_time_diff = Int(25 - min_down_time)
                for j in 1:down_time_diff
                    z_range = j:(Int(min_down_time + j - 1))
                    subgraph_index = 26 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (down_time_diff + 1):25
                    z_range = j:25
                    subgraph_index = 26 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 2 # s_list is length 28
            if min_down_time > 24
                down_time_diff = Int(49 - min_down_time)

                for j in 1:down_time_diff
                    z_range = j:(Int(min_down_time + j - 1))
                    subgraph_index = 26 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(DA_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (down_time_diff + 1):25
                    z_range = j:49
                    subgraph_index = 26 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(DA_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:25
                    subgraph_index = 26 - j
                    z_range = j:(Int(min_down_time + j - 1))
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(DA_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        else
            for j in 1:25
                subgraph_index = 26 - j
                z_range = j:(Int(min_down_time + j - 1))
                @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(DA_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
            end
        end
    end

    num_DA_gen = size(gen_DA, 1)
    for j in 1:num_DA_gen
        gen_name = gen_DA[j, "Generator Name"]
        bus_name = gen_DA[j, "Bus of Connection"]

        min_up_time = gen_DA[j, "Min Up Time (h)"]
        min_down_time = gen_DA[j, "Min Down Time (h)"]

        max_ramp_up = gen_DA[j, "Max Ramp Up (MW/min)"] * 60 * dt
        max_ramp_down = gen_DA[j, "Max Ramp Down (MW/min)"] * 60 * dt

        max_cap = gen_DA[j, "Max Capacity (MW)"]
        min_cap = gen_DA[j, "Min Stable Level (MW)"]

        max_startup = gen_DA[j, "Max Capacity (MW)"]
        max_shutdown = gen_DA[j, "Max Capacity (MW)"]

        s_list = get_s_list(DA_graph_set, bus_name, gen_name, day_num)
        z_list = get_z_list(DA_graph_set, bus_name, gen_name, day_num)
        if !(isequal(0, min_up_time))
            link_subgraphs_uptime(DA_graph, DA_graph_set, min_up_time, bus_name, gen_name, s_list, day_num)
        end
        if !(isequal(0, min_down_time))
            link_subgraphs_downtime(DA_graph, DA_graph_set, min_down_time, bus_name, gen_name, z_list, day_num)
        end

        if day_num != 1
            if day_num == 2
                next_DA_graph_index = 2
                last_DA_graph_index = 1
            else
                next_DA_graph_index = 3
                last_DA_graph_index = 2
            end
            bus_next = getsubgraphs(DA_graph_set[next_DA_graph_index])[1][:bus][bus_name]
            bus_last = getsubgraphs(DA_graph_set[last_DA_graph_index])[24][:bus][bus_name]
            bus_last_HA = getsubgraphs(HA_graph)[1][:bus][bus_name]
            @linkconstraint(DA_graph, bus_next[:x][gen_name] - value(bus_last[:x][gen_name]) == bus_next[:s][gen_name] - bus_next[:z][gen_name])

            @linkconstraint(DA_graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                value(bus_last_HA[:Gen_p][gen_name]) - value(bus_last_HA[:Gen_m][gen_name]) <=
                (max_startup - max_ramp_up - min_cap) * bus_next[:s][gen_name] +
                (max_ramp_up + min_cap) * bus_next[:x][gen_name] - min_cap * value(bus_last[:x][gen_name])
            )

            @linkconstraint(DA_graph, value(bus_last_HA[:Gen_p][gen_name]) + value(bus_last_HA[:Gen_m][gen_name]) -
                bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                (max_shutdown - max_ramp_down - min_cap) * bus_next[:z][gen_name] +
                (max_ramp_down + min_cap) * value(bus_last[:x][gen_name]) - min_cap * bus_next[:x][gen_name]
            )
        end
    end
end

function link_sol_over_HA_new_day(graph, graph_last; dt = .25, monolithic = true)

    #This function needs to link the new day to the HA solution of the last day

    HA_subgraph_next = getsubgraphs(graph)[10]
    HA_subgraph_last = getsubgraphs(graph_last)[105]
    ST_subgraph_next = getsubgraphs(graph)[2]
    ST_subgraph_last = getsubgraphs(graph_last)[9]
    DA_subgraph_next = getsubgraphs(graph)[1]
    DA_subgraph_last = getsubgraphs(graph_last)[1]

    num_conv_gen = size(gen_data_conv, 1)

    DA_len = size(gen_DA, 1)


    for k in 1:DA_len
        gen_name = gen_DA[k, "Generator Name"]
        bus_name = gen_DA[k, "Bus of Connection"]

        gen_max_cap = gen_DA[k, "Max Capacity (MW)"]
        gen_min_cap = gen_DA[k, "Min Stable Level (MW)"]

        max_ramp_up = gen_DA[k, "Max Ramp Up (MW/min)"] * 60 * dt# * dt #TODO: check if this should be multiplied by dt or not
        max_ramp_down = gen_DA[k, "Max Ramp Down (MW/min)"] * 60 * dt# * dt

        max_startup = gen_DA[k, "Max Capacity (MW)"]
        max_shutdown = gen_DA[k, "Max Capacity (MW)"]

        bus_next = getsubgraphs(HA_subgraph_next)[1][:bus][bus_name]
        bus_last = getsubgraphs(HA_subgraph_last)[1][:bus][bus_name]
        bus_last_gen_p = value(bus_last[:Gen_p][gen_name])
        bus_last_gen_m = value(bus_last[:Gen_m][gen_name])


        bus_next_DA = getsubgraphs(DA_subgraph_next)[1][:bus][bus_name]
        bus_last_DA = getsubgraphs(DA_subgraph_last)[24][:bus][bus_name]

        if monolithic
            DA_x_next = bus_next_DA[:x][gen_name]
            DA_s_next = bus_next_DA[:s][gen_name]
            DA_z_next = bus_next_DA[:z][gen_name]
            DA_x_last = value(bus_last_DA[:x][gen_name])
        else
            DA_x_next = value(bus_next_DA[:x][gen_name])
            DA_s_next = value(bus_next_DA[:s][gen_name])
            DA_z_next = value(bus_next_DA[:z][gen_name])
            DA_x_last = value(bus_last_DA[:x][gen_name])
        end


        @linkconstraint(HA_subgraph_next, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
            bus_last_gen_p - bus_last_gen_m <=
            (max_startup - max_ramp_up - gen_min_cap) * DA_s_next +
            (max_ramp_up + gen_min_cap) * DA_x_next - gen_min_cap * DA_x_last
        )

        @linkconstraint(HA_subgraph_next, bus_last_gen_p + bus_last_gen_m -
            bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
            (max_shutdown - max_ramp_down - gen_min_cap) * DA_z_next +
            (max_ramp_down + gen_min_cap) * DA_x_last - gen_min_cap * DA_x_next
        )
    end

    ST_len = size(gen_ST, 1)
    for k in 1:ST_len
        gen_name = gen_ST[k, "Generator Name"]
        bus_name = gen_ST[k, "Bus of Connection"]

        gen_max_cap = gen_ST[k, "Max Capacity (MW)"]
        gen_min_cap = gen_ST[k, "Min Stable Level (MW)"]

        max_ramp_up = gen_ST[k, "Max Ramp Up (MW/min)"] * 60 * dt
        max_ramp_down = gen_ST[k, "Max Ramp Down (MW/min)"] * 60 * dt

        max_startup = gen_ST[k, "Max Capacity (MW)"]
        max_shutdown = gen_ST[k, "Max Capacity (MW)"]

        bus_next = getsubgraphs(HA_subgraph_next)[1][:bus][bus_name]
        bus_last = getsubgraphs(HA_subgraph_last)[1][:bus][bus_name]
        bus_last_gen_p = value(bus_last[:Gen_p][gen_name])
        bus_last_gen_m = value(bus_last[:Gen_m][gen_name])

        bus_next_ST = getsubgraphs(ST_subgraph_next)[1][:bus][bus_name]
        bus_last_ST = getsubgraphs(ST_subgraph_last)[12][:bus][bus_name]

        if monolithic
            ST_x_next = bus_next_ST[:x][gen_name]
            ST_s_next = bus_next_ST[:s][gen_name]
            ST_z_next = bus_next_ST[:z][gen_name]
            ST_x_last = value(bus_last_ST[:x][gen_name])
        else
            ST_x_next = value(bus_next_ST[:x][gen_name])
            ST_s_next = value(bus_next_ST[:s][gen_name])
            ST_z_next = value(bus_next_ST[:z][gen_name])
            ST_x_last = value(bus_last_ST[:x][gen_name])
        end

        @linkconstraint(HA_subgraph_next, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
            bus_last_gen_p - bus_last_gen_m <=
            (max_startup - max_ramp_up - gen_min_cap) * ST_s_next +
            (max_ramp_up + gen_min_cap) * ST_x_next - gen_min_cap * ST_x_last
        )

        @linkconstraint(HA_subgraph_next, bus_last_gen_p + bus_last_gen_m -
            bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
            (max_shutdown - max_ramp_down - gen_min_cap) * ST_z_next +
            (max_ramp_down + gen_min_cap) * ST_x_last - gen_min_cap * ST_x_next
        )
    end
end
