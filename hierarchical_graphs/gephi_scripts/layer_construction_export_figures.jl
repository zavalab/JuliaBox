function link_between_ST_no_sol(graph; dt = .25, day_num = 1, graph_last = nothing)

    function get_s_list(ST_graphs, bus_name, gen_name, i)
        if i == 1
            s_list = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1]
        elseif i == 2
            s_list1 = [getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1]
            s_list2 = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list = vcat(s_list1, s_list2)
        elseif i == 3
            s_list1 = [getsubgraphs(ST_graphs[3])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1]
            s_list2 = [getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list3 = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list = vcat(s_list1, s_list2, s_list3)
        else
            s_list1 = [getsubgraphs(ST_graphs[i])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1]
            s_list2 = [getsubgraphs(ST_graphs[i - 1])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list3 = [getsubgraphs(ST_graphs[i - 2])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list4 = [getsubgraphs(ST_graphs[i - 3])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list = vcat(s_list1, s_list2, s_list3, s_list4)
        end
        return s_list
    end

    function get_z_list(ST_graphs, bus_name, gen_name, i)
        if i == 1
            z_list = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1]
        elseif i == 2
            z_list1 = [getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1]
            z_list2 = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list = vcat(z_list1, z_list2)
        elseif i == 3
            z_list1 = [getsubgraphs(ST_graphs[3])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1]
            z_list2 = [getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list3 = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list = vcat(z_list1, z_list2, z_list3)
        else
            z_list1 = [getsubgraphs(ST_graphs[i])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1]
            z_list2 = [getsubgraphs(ST_graphs[i - 1])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list3 = [getsubgraphs(ST_graphs[i - 2])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list4 = [getsubgraphs(ST_graphs[i - 3])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list = vcat(z_list1, z_list2, z_list3, z_list4)
        end
        return z_list
    end

    function link_subgraphs_uptime(graph, ST_graphs, min_up_time, bus_name, gen_name, s_list, i)
        if i == 1
            if min_up_time >= 4
                for j in 1:15
                    subgraph_index = 17 - j
                    s_range = (j + 1):16
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                #@linkconstraint(graph, sum(s_list[2:16]) <= getsubgraphs(ST_graphs)[16][:bus][bus_name][:x][gen_name])
            else
                up_time_diff = Int(4 * (4 - min_up_time))
                for j in 1:up_time_diff
                    s_range = (j + 1):(Int(min_up_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):15
                    s_range = (j + 1):16
                    subgraph_index = 17 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 2 # s_list is length 28
            if (min_up_time >= 3) && (min_up_time < 7)
                up_time_diff = Int(4 * (7 - min_up_time))

                for j in 1:up_time_diff
                    s_range = (j + 1):(Int(min_up_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):16
                    s_range = (j + 1):28
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            elseif (min_up_time >= 7)
                for j in 1:16
                    subgraph_index = 17 - j
                    s_range = (j + 1):28
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:16
                    subgraph_index = 17 - j
                    s_range = (j + 1):(Int(min_up_time * 4 + j - 1))
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 3 #s_list is length 40
            if (min_up_time >= 6)
                up_time_diff = Int(4 * (10 - min_up_time))
                for j in 1:up_time_diff
                    s_range = (j + 1):(Int(min_up_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):16
                    s_range = (j + 1):40
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:16
                    subgraph_index = 17 - j
                    s_range = (j + 1):(Int(min_up_time * 4 + j - 1))
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        else #s_list is length 52
            for j in 1:16
                subgraph_index = 17 - j
                s_range = (j + 1):(Int(min_up_time * 4 + j - 1))
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
                #@linkconstraint(graph, sum(z_list[1:16]) <= 1 - getsubgraphs(ST_graphs)[16][:bus][bus_name][:x][gen_name])
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
                    @linkconstraint(graph, sum(z_list[z_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 3 #s_list is length 40
            if (min_down_time >= 6)
                down_time_diff = Int(4 * (10 - min_down_time))
                for j in 1:down_time_diff
                    z_range = j:(Int(min_down_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (down_time_diff + 1):16
                    z_range = j:40
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:16
                    subgraph_index = 17 - j
                    z_range = j:(Int(min_down_time * 4 + j - 1))
                    @linkconstraint(graph, sum(z_list[z_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        else #s_list is length 52
            for j in 1:16
                subgraph_index = 17 - j
                z_range = j:(Int(min_down_time * 4 + j - 1))
                @linkconstraint(graph, sum(z_list[z_range]) <= getsubgraphs(ST_graphs[i])[subgraph_index][:bus][bus_name][:x][gen_name])
            end
        end
    end

    num_ST_gen = size(gen_ST, 1)

    for i in 1:8
        if day_num != 1
            if i == 1
                ST_subgraph_set1 = getsubgraphs(graph_last)[7:9]
                ST_subgraph_set2 = getsubgraphs(graph)[2:2]
                ST_graphs = vcat(ST_subgraph_set1, ST_subgraph_set2)
            elseif i == 2
                ST_subgraph_set1 = getsubgraphs(graph_last)[8:9]
                ST_subgraph_set2 = getsubgraphs(graph)[2:3]
                ST_graphs = vcat(ST_subgraph_set1, ST_subgraph_set2)
            elseif i == 3
                ST_subgraph_set1 = getsubgraphs(graph_last)[9:9]
                ST_subgraph_set2 = getsubgraphs(graph)[2:4]
                ST_graphs = vcat(ST_subgraph_set1, ST_subgraph_set2)
            else
                ST_graphs = getsubgraphs(graph)[2:(i + 1)]
            end
        else
            ST_graphs = getsubgraphs(graph)[2:9]
        end
        ST_graphs_len = length(ST_graphs)

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
                link_subgraphs_uptime(graph, ST_graphs, min_up_time, bus_name, gen_name, s_list, ST_graphs_len)
            end
            if !(isequal(0, min_down_time))
                link_subgraphs_downtime(graph, ST_graphs, min_down_time, bus_name, gen_name, z_list, ST_graphs_len)
            end

            ST_graphs = getsubgraphs(graph)[2:9]

            if i != 1
                bus_next = getsubgraphs(ST_graphs[i])[1][:bus][bus_name]
                bus_last = getsubgraphs(ST_graphs[i - 1])[12][:bus][bus_name]
                HA_subgraph = getsubgraphs(graph)[9 + (i - 1) * 12]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]
                @linkconstraint(graph, bus_next[:x][gen_name] - bus_last[:x][gen_name] == bus_next[:s][gen_name] - bus_next[:z][gen_name])

                @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                    bus_last_HA[:Gen_p][gen_name] - bus_last_HA[:Gen_m][gen_name] <=
                    (max_startup - max_ramp_up - min_cap) * bus_next[:s][gen_name] +
                    (max_ramp_up + min_cap) * bus_next[:x][gen_name] - min_cap * bus_last[:x][gen_name]
                )

                @linkconstraint(graph, bus_last_HA[:Gen_p][gen_name] + bus_last_HA[:Gen_m][gen_name] -
                    bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                    (max_shutdown - max_ramp_down - min_cap) * bus_next[:z][gen_name] +
                    (max_ramp_down + min_cap) * bus_last[:x][gen_name] - min_cap * bus_next[:x][gen_name]
                )
            elseif day_num != 1
                bus_next = getsubgraphs(ST_graphs[1])[1][:bus][bus_name]
                bus_last = getsubgraphs(getsubgraphs(graph_last)[9])[12][:bus][bus_name]
                HA_subgraph = getsubgraphs(graph_last)[105]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]
                @linkconstraint(graph, bus_next[:x][gen_name] - bus_last[:x][gen_name] == bus_next[:s][gen_name] - bus_next[:z][gen_name])

                @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                    bus_last_HA[:Gen_p][gen_name] - bus_last_HA[:Gen_m][gen_name] <=
                    max_startup - max_ramp_up - min_cap * bus_next[:s][gen_name] +
                    (max_ramp_up + min_cap) * bus_next[:x][gen_name] - min_cap * bus_last[:x][gen_name]
                )

                @linkconstraint(graph, bus_last_HA[:Gen_p][gen_name] + bus_last_HA[:Gen_m][gen_name] -
                    bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                    (max_shutdown - max_ramp_down - min_cap) * bus_next[:z][gen_name] +
                    (max_ramp_down + min_cap) * bus_last[:x][gen_name] - min_cap * bus_next[:x][gen_name]
                )
            end
        end
    end

    num_DA_gen = size(gen_DA, 1)
    ST_graphs = getsubgraphs(graph)[2:9]

    if day_num == 1
        DA_iter_index = 2:8
    else
        DA_iter_index = 1:8
    end

    for i in DA_iter_index
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
                bus_next = getsubgraphs(ST_graphs[i])[1][:bus][bus_name]
                bus_last = getsubgraphs(getsubgraphs(graph_last)[9])[12][:bus][bus_name]

                HA_subgraph = getsubgraphs(graph_last)[105]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]

                DA_subgraph_last = getsubgraphs(graph_last)[1]
                bus_last_DA = getsubgraphs(DA_subgraph_last)[24][:bus][bus_name]

                DA_subgraph_next = getsubgraphs(graph)[1]
                bus_next_DA = getsubgraphs(DA_subgraph_next)[1][:bus][bus_name]
            else
                bus_next = getsubgraphs(ST_graphs[i])[1][:bus][bus_name]
                bus_last = getsubgraphs(ST_graphs[i - 1])[12][:bus][bus_name]

                HA_subgraph = getsubgraphs(graph)[9 + (i - 1) * 12]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]

                DA_subgraph = getsubgraphs(graph)[1]
                bus_last_DA = getsubgraphs(DA_subgraph)[(i - 1) * 3][:bus][bus_name]

                DA_subgraph = getsubgraphs(graph)[1]
                bus_next_DA = getsubgraphs(DA_subgraph)[1 + (i - 1) * 3][:bus][bus_name]

            end

            @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                bus_last_HA[:Gen_p][gen_name] - bus_last_HA[:Gen_m][gen_name] <=
                (max_startup - max_ramp_up / dt - min_cap) * bus_next_DA[:s][gen_name] +
                (max_ramp_up / dt + min_cap) * bus_next_DA[:x][gen_name] - min_cap * bus_last_DA[:x][gen_name]
            )

            @linkconstraint(graph, bus_last_HA[:Gen_p][gen_name] + bus_last_HA[:Gen_m][gen_name] -
                bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                (max_shutdown - max_ramp_down / dt - min_cap) * bus_next_DA[:z][gen_name] +
                (max_ramp_down / dt + min_cap) * bus_last_DA[:x][gen_name] - min_cap * bus_next_DA[:x][gen_name]
            )
        end
    end
end

function link_between_DA_no_sol(DA_graph, DA_graph_set, day_num; dt = 1, HA_graph = nothing)

    function get_s_list(DA_graph_set, bus_name, gen_name, day_num)

        if day_num == 1
            s_list = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:s][gen_name] for j in 25:-1:1]
        elseif day_num == 2
            s_list1 = [getsubgraphs(DA_graph_set[2])[j][:bus][bus_name][:s][gen_name] for j in 25:-1:1]
            s_list2 = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:s][gen_name] for j in 24:-1:1]
            s_list = vcat(s_list1, s_list2)
        else
            s_list1 = [getsubgraphs(DA_graph_set[3])[j][:bus][bus_name][:s][gen_name] for j in 25:-1:1]
            s_list2 = [getsubgraphs(DA_graph_set[2])[j][:bus][bus_name][:s][gen_name] for j in 24:-1:1]
            s_list3 = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:s][gen_name] for j in 24:-1:1]
            s_list = vcat(s_list1, s_list2, s_list3)
        end
        return s_list
    end

    function get_z_list(DA_graph_set, bus_name, gen_name, day_num)
        if day_num == 1
            z_list = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:z][gen_name] for j in 25:-1:1]
        elseif day_num == 2
            z_list1 = [getsubgraphs(DA_graph_set[2])[j][:bus][bus_name][:z][gen_name] for j in 25:-1:1]
            z_list2 = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:z][gen_name] for j in 24:-1:1]
            z_list = vcat(z_list1, z_list2)
        else
            z_list1 = [getsubgraphs(DA_graph_set[3])[j][:bus][bus_name][:z][gen_name] for j in 25:-1:1]
            z_list2 = [getsubgraphs(DA_graph_set[2])[j][:bus][bus_name][:z][gen_name] for j in 24:-1:1]
            z_list3 = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:z][gen_name] for j in 24:-1:1]
            z_list = vcat(z_list1, z_list2, z_list3)
        end
        return z_list
    end

    function link_subgraphs_uptime(graph, DA_graphs, min_up_time, bus_name, gen_name, s_list, i)
        if i == 1# Look at code for elseif i == 3 below
            if min_up_time >24
                for j in 1:24
                    subgraph_index = 26 - j
                    s_range = (j + 1):25
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                #@linkconstraint(graph, sum(s_list[2:16]) <= getsubgraphs(ST_graphs)[16][:bus][bus_name][:x][gen_name])
            else
                up_time_diff = Int(25 - min_up_time)
                for j in 1:up_time_diff
                    s_range = (j + 1):(Int(min_up_time + j - 1))
                    subgraph_index = 26 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):24
                    s_range = (j + 1):25
                    subgraph_index = 26 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 2 # s_list is length 28
            if min_up_time > 24
                up_time_diff = Int(49 - min_up_time)

                for j in 1:up_time_diff
                    s_range = (j + 1):(Int(min_up_time + j - 1))
                    subgraph_index = 26 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):25
                    s_range = (j + 1):49
                    subgraph_index = 26 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:25
                    subgraph_index = 26 - j
                    s_range = (j + 1):(Int(min_up_time + j - 1))
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        else
            for j in 1:25
                subgraph_index = 26 - j
                s_range = (j + 1):(Int(min_up_time + j - 1))
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
                @linkconstraint(graph, sum(z_list[z_range]) <= getsubgraphs(DA_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
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
            @linkconstraint(DA_graph, bus_next[:x][gen_name] - bus_last[:x][gen_name] == bus_next[:s][gen_name] - bus_next[:z][gen_name])

            @linkconstraint(DA_graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                bus_last_HA[:Gen_p][gen_name] - bus_last_HA[:Gen_m][gen_name] <=
                (max_startup - max_ramp_up - min_cap) * bus_next[:s][gen_name] +
                (max_ramp_up + min_cap) * bus_next[:x][gen_name] - min_cap * bus_last[:x][gen_name]
            )

            @linkconstraint(DA_graph, bus_last_HA[:Gen_p][gen_name] + bus_last_HA[:Gen_m][gen_name] -
                bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                (max_shutdown - max_ramp_down - min_cap) * bus_next[:z][gen_name] +
                (max_ramp_down + min_cap) * bus_last[:x][gen_name] - min_cap * bus_next[:x][gen_name]
            )
        end
    end
end

function link_sol_over_HA_new_day_no_sol(graph, graph_last; dt = .25, monolithic = true)

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
        bus_last_gen_p = bus_last[:Gen_p][gen_name]
        bus_last_gen_m = bus_last[:Gen_m][gen_name]


        bus_next_DA = getsubgraphs(DA_subgraph_next)[1][:bus][bus_name]
        bus_last_DA = getsubgraphs(DA_subgraph_last)[24][:bus][bus_name]

        if monolithic
            DA_x_next = bus_next_DA[:x][gen_name]
            DA_s_next = bus_next_DA[:s][gen_name]
            DA_z_next = bus_next_DA[:z][gen_name]
            DA_x_last = bus_last_DA[:x][gen_name]
        else
            DA_x_next = bus_next_DA[:x][gen_name]
            DA_s_next = bus_next_DA[:s][gen_name]
            DA_z_next = bus_next_DA[:z][gen_name]
            DA_x_last = bus_last_DA[:x][gen_name]
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
        bus_last_gen_p = bus_last[:Gen_p][gen_name]
        bus_last_gen_m = bus_last[:Gen_m][gen_name]

        bus_next_ST = getsubgraphs(ST_subgraph_next)[1][:bus][bus_name]
        bus_last_ST = getsubgraphs(ST_subgraph_last)[12][:bus][bus_name]

        if monolithic
            ST_x_next = bus_next_ST[:x][gen_name]
            ST_s_next = bus_next_ST[:s][gen_name]
            ST_z_next = bus_next_ST[:z][gen_name]
            ST_x_last = bus_last_ST[:x][gen_name]
        else
            ST_x_next = bus_next_ST[:x][gen_name]
            ST_s_next = bus_next_ST[:s][gen_name]
            ST_z_next = bus_next_ST[:z][gen_name]
            ST_x_last = bus_last_ST[:x][gen_name]
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


function link_between_ST_4hr(graph; dt = .25, day_num = 1, graph_last = nothing)

    function get_s_list(ST_graphs, bus_name, gen_name, i)
        if i == 1
            s_list = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1]
        elseif i == 2
            s_list1 = [getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1]
            s_list2 = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list = vcat(s_list1, s_list2)
        elseif i == 3
            s_list1 = [getsubgraphs(ST_graphs[3])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1]
            s_list2 = [getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list3 = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list = vcat(s_list1, s_list2, s_list3)
        else
            s_list1 = [getsubgraphs(ST_graphs[i])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1]
            s_list2 = [getsubgraphs(ST_graphs[i - 1])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list3 = [getsubgraphs(ST_graphs[i - 2])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list4 = [getsubgraphs(ST_graphs[i - 3])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list = vcat(s_list1, s_list2, s_list3, s_list4)
        end
        return s_list
    end

    function get_z_list(ST_graphs, bus_name, gen_name, i)
        if i == 1
            z_list = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1]
        elseif i == 2
            z_list1 = [getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1]
            z_list2 = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list = vcat(z_list1, z_list2)
        elseif i == 3
            z_list1 = [getsubgraphs(ST_graphs[3])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1]
            z_list2 = [getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list3 = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list = vcat(z_list1, z_list2, z_list3)
        else
            z_list1 = [getsubgraphs(ST_graphs[i])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1]
            z_list2 = [getsubgraphs(ST_graphs[i - 1])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list3 = [getsubgraphs(ST_graphs[i - 2])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list4 = [getsubgraphs(ST_graphs[i - 3])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list = vcat(z_list1, z_list2, z_list3, z_list4)
        end
        return z_list
    end

    function link_subgraphs_uptime(graph, ST_graphs, min_up_time, bus_name, gen_name, s_list, i)
        if i == 1
            if min_up_time >= 4
                for j in 1:15
                    subgraph_index = 17 - j
                    s_range = (j + 1):16
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                #@linkconstraint(graph, sum(s_list[2:16]) <= getsubgraphs(ST_graphs)[16][:bus][bus_name][:x][gen_name])
            else
                up_time_diff = Int(4 * (4 - min_up_time))
                for j in 1:up_time_diff
                    s_range = (j + 1):(Int(min_up_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):15
                    s_range = (j + 1):16
                    subgraph_index = 17 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 2 # s_list is length 28
            if (min_up_time >= 3) && (min_up_time < 7)
                up_time_diff = Int(4 * (7 - min_up_time))

                for j in 1:up_time_diff
                    s_range = (j + 1):(Int(min_up_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):16
                    s_range = (j + 1):28
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            elseif (min_up_time >= 7)
                for j in 1:16
                    subgraph_index = 17 - j
                    s_range = (j + 1):28
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:16
                    subgraph_index = 17 - j
                    s_range = (j + 1):(Int(min_up_time * 4 + j - 1))
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 3 #s_list is length 40
            if (min_up_time >= 6)
                up_time_diff = Int(4 * (10 - min_up_time))
                for j in 1:up_time_diff
                    s_range = (j + 1):(Int(min_up_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):16
                    s_range = (j + 1):40
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:16
                    subgraph_index = 17 - j
                    s_range = (j + 1):(Int(min_up_time * 4 + j - 1))
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        else #s_list is length 52
            for j in 1:16
                subgraph_index = 17 - j
                s_range = (j + 1):(Int(min_up_time * 4 + j - 1))
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
                #@linkconstraint(graph, sum(z_list[1:16]) <= 1 - getsubgraphs(ST_graphs)[16][:bus][bus_name][:x][gen_name])
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
                    @linkconstraint(graph, sum(z_list[z_range]) <= getsubgraphs(ST_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 3 #s_list is length 40
            if (min_down_time >= 6)
                down_time_diff = Int(4 * (10 - min_down_time))
                for j in 1:down_time_diff
                    z_range = j:(Int(min_down_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (down_time_diff + 1):16
                    z_range = j:40
                    subgraph_index = 17 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:16
                    subgraph_index = 17 - j
                    z_range = j:(Int(min_down_time * 4 + j - 1))
                    @linkconstraint(graph, sum(z_list[z_range]) <= getsubgraphs(ST_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        else #s_list is length 52
            for j in 1:16
                subgraph_index = 17 - j
                z_range = j:(Int(min_down_time * 4 + j - 1))
                @linkconstraint(graph, sum(z_list[z_range]) <= getsubgraphs(ST_graphs[i])[subgraph_index][:bus][bus_name][:x][gen_name])
            end
        end
    end

    num_ST_gen = size(gen_ST, 1)

    for i in 1:1
        if day_num != 1
            if i == 1
                ST_subgraph_set1 = getsubgraphs(graph_last)[7:9]
                ST_subgraph_set2 = getsubgraphs(graph)[2:2]
                ST_graphs = vcat(ST_subgraph_set1, ST_subgraph_set2)
            elseif i == 2
                ST_subgraph_set1 = getsubgraphs(graph_last)[8:9]
                ST_subgraph_set2 = getsubgraphs(graph)[2:3]
                ST_graphs = vcat(ST_subgraph_set1, ST_subgraph_set2)
            elseif i == 3
                ST_subgraph_set1 = getsubgraphs(graph_last)[9:9]
                ST_subgraph_set2 = getsubgraphs(graph)[2:4]
                ST_graphs = vcat(ST_subgraph_set1, ST_subgraph_set2)
            else
                ST_graphs = getsubgraphs(graph)[2:(i + 1)]
            end
        else
            ST_graphs = getsubgraphs(graph)[2:9]
        end
        ST_graphs_len = length(ST_graphs)

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
                link_subgraphs_uptime(graph, ST_graphs, min_up_time, bus_name, gen_name, s_list, ST_graphs_len)
            end
            if !(isequal(0, min_down_time))
                link_subgraphs_downtime(graph, ST_graphs, min_down_time, bus_name, gen_name, z_list, ST_graphs_len)
            end

            ST_graphs = getsubgraphs(graph)[2:9]

            if i != 1
                bus_next = getsubgraphs(ST_graphs[i])[1][:bus][bus_name]
                bus_last = getsubgraphs(ST_graphs[i - 1])[12][:bus][bus_name]
                HA_subgraph = getsubgraphs(graph)[9 + (i - 1) * 12]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]
                @linkconstraint(graph, bus_next[:x][gen_name] - bus_last[:x][gen_name] == bus_next[:s][gen_name] - bus_next[:z][gen_name])

                @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                    bus_last_HA[:Gen_p][gen_name] - bus_last_HA[:Gen_m][gen_name] <=
                    (max_startup - max_ramp_up - min_cap) * bus_next[:s][gen_name] +
                    (max_ramp_up + min_cap) * bus_next[:x][gen_name] - min_cap * bus_last[:x][gen_name]
                )

                @linkconstraint(graph, bus_last_HA[:Gen_p][gen_name] + bus_last_HA[:Gen_m][gen_name] -
                    bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                    (max_shutdown - max_ramp_down - min_cap) * bus_next[:z][gen_name] +
                    (max_ramp_down + min_cap) * bus_last[:x][gen_name] - min_cap * bus_next[:x][gen_name]
                )
            elseif day_num != 1
                bus_next = getsubgraphs(ST_graphs[1])[1][:bus][bus_name]
                bus_last = getsubgraphs(getsubgraphs(graph_last)[9])[12][:bus][bus_name]
                HA_subgraph = getsubgraphs(graph_last)[105]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]
                @linkconstraint(graph, bus_next[:x][gen_name] - bus_last[:x][gen_name] == bus_next[:s][gen_name] - bus_next[:z][gen_name])

                @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                    bus_last_HA[:Gen_p][gen_name] - bus_last_HA[:Gen_m][gen_name] <=
                    max_startup - max_ramp_up - min_cap * bus_next[:s][gen_name] +
                    (max_ramp_up + min_cap) * bus_next[:x][gen_name] - min_cap * bus_last[:x][gen_name]
                )

                @linkconstraint(graph, bus_last_HA[:Gen_p][gen_name] + bus_last_HA[:Gen_m][gen_name] -
                    bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                    (max_shutdown - max_ramp_down - min_cap) * bus_next[:z][gen_name] +
                    (max_ramp_down + min_cap) * bus_last[:x][gen_name] - min_cap * bus_next[:x][gen_name]
                )
            end
        end
    end
end


function link_between_DA_4hr(DA_graph, DA_graph_set, day_num; dt = 1, HA_graph = nothing)

    function get_s_list(DA_graph_set, bus_name, gen_name, day_num)

        if day_num == 1
            s_list = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:s][gen_name] for j in 4:-1:1]
        elseif day_num == 2
            s_list1 = [getsubgraphs(DA_graph_set[2])[j][:bus][bus_name][:s][gen_name] for j in 25:-1:1]
            s_list2 = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:s][gen_name] for j in 24:-1:1]
            s_list = vcat(s_list1, s_list2)
        else
            s_list1 = [getsubgraphs(DA_graph_set[3])[j][:bus][bus_name][:s][gen_name] for j in 25:-1:1]
            s_list2 = [getsubgraphs(DA_graph_set[2])[j][:bus][bus_name][:s][gen_name] for j in 24:-1:1]
            s_list3 = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:s][gen_name] for j in 24:-1:1]
            s_list = vcat(s_list1, s_list2, s_list3)
        end
        return s_list
    end

    function get_z_list(DA_graph_set, bus_name, gen_name, day_num)
        if day_num == 1
            z_list = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:z][gen_name] for j in 25:-1:1]
        elseif day_num == 2
            z_list1 = [getsubgraphs(DA_graph_set[2])[j][:bus][bus_name][:z][gen_name] for j in 25:-1:1]
            z_list2 = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:z][gen_name] for j in 24:-1:1]
            z_list = vcat(z_list1, z_list2)
        else
            z_list1 = [getsubgraphs(DA_graph_set[3])[j][:bus][bus_name][:z][gen_name] for j in 25:-1:1]
            z_list2 = [getsubgraphs(DA_graph_set[2])[j][:bus][bus_name][:z][gen_name] for j in 24:-1:1]
            z_list3 = [getsubgraphs(DA_graph_set[1])[j][:bus][bus_name][:z][gen_name] for j in 24:-1:1]
            z_list = vcat(z_list1, z_list2, z_list3)
        end
        return z_list
    end

    function link_subgraphs_uptime(graph, DA_graphs, min_up_time, bus_name, gen_name, s_list, i)
        if i == 1# Look at code for elseif i == 3 below
            if min_up_time >3
                for j in 1:3
                    subgraph_index = 5 - j
                    s_range = (j + 1):4
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                #@linkconstraint(graph, sum(s_list[2:16]) <= getsubgraphs(ST_graphs)[16][:bus][bus_name][:x][gen_name])
            else
                up_time_diff = Int(3 - min_up_time)
                for j in 1:up_time_diff
                    s_range = (j + 1):(Int(min_up_time + j - 1))
                    subgraph_index = 5 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):3
                    s_range = (j + 1):4
                    subgraph_index = 5 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 2 # s_list is length 28
            if min_up_time > 24
                up_time_diff = Int(49 - min_up_time)

                for j in 1:up_time_diff
                    s_range = (j + 1):(Int(min_up_time + j - 1))
                    subgraph_index = 26 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):25
                    s_range = (j + 1):49
                    subgraph_index = 26 - j
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                for j in 1:25
                    subgraph_index = 26 - j
                    s_range = (j + 1):(Int(min_up_time + j - 1))
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[2])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        else
            for j in 1:25
                subgraph_index = 26 - j
                s_range = (j + 1):(Int(min_up_time + j - 1))
                @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
            end
        end
    end

    function link_subgraphs_downtime(graph, DA_graphs, min_down_time, bus_name, gen_name, z_list, i)
        if i == 1
            if min_down_time > 3
                for j in 1:3
                    subgraph_index = 5 - j
                    z_range = j:4
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                down_time_diff = Int(4 - min_down_time)
                for j in 1:down_time_diff
                    z_range = j:(Int(min_down_time + j - 1))
                    subgraph_index =5 - j
                    @linkconstraint(graph, sum(z_list[z_range]) <= 1 - getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (down_time_diff + 1):4
                    z_range = j:4
                    subgraph_index = 5 - j
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
                @linkconstraint(graph, sum(z_list[z_range]) <= getsubgraphs(DA_graphs[3])[subgraph_index][:bus][bus_name][:x][gen_name])
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
            @linkconstraint(DA_graph, bus_next[:x][gen_name] - bus_last[:x][gen_name] == bus_next[:s][gen_name] - bus_next[:z][gen_name])

            @linkconstraint(DA_graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                bus_last_HA[:Gen_p][gen_name] - bus_last_HA[:Gen_m][gen_name] <=
                (max_startup - max_ramp_up - min_cap) * bus_next[:s][gen_name] +
                (max_ramp_up + min_cap) * bus_next[:x][gen_name] - min_cap * bus_last[:x][gen_name]
            )

            @linkconstraint(DA_graph, bus_last_HA[:Gen_p][gen_name] + bus_last_HA[:Gen_m][gen_name] -
                bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                (max_shutdown - max_ramp_down - min_cap) * bus_next[:z][gen_name] +
                (max_ramp_down + min_cap) * bus_last[:x][gen_name] - min_cap * bus_next[:x][gen_name]
            )
        end
    end
end
