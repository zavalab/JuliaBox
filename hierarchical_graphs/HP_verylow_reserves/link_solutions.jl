using Plasmo
using DelimitedFiles, CSV, DataFrames
using Gurobi

"""
    link_DA_sol_to_ST(graph::OptiGraph, ST_subgraph:OptiGraph, tau_0, dt, gen_com_set)

Links the solutions in the DA level to the ST level for a given ST level subgraph.
 * graph - full day OptiGraph with all three layers (should have DAUC as first Optigraph)
 * ST_subgraph - ST_subgraph being linked to DA layer
 * tau_0 - set of hour time points corresponding to the ST graph (1-4, 4-7, 7-10, etc.)
 * dt - discretization size (in hours)
 * gen_com_set - set of generators to be committed
"""
function link_DA_sol_to_ST(graph::OptiGraph, ST_subgraph, tau_0 = 1:4, dt = .25, gen_com_set = gen_ST)

    # Get the DAUC subgraph
    DA_subgraph = getsubgraphs(graph)[1]
    # get the number of DA generators
    num_DA_gen = size(gen_DA, 1)

    # Iterate through the DA-UC time points that correspond to the problem
    for (i, index) in enumerate(tau_0)
        if index <= length(getsubgraphs(DA_subgraph))

            # Get the tiem point subgraph from DA_subgraph
            DA_subgraph_i = getsubgraphs(DA_subgraph)[index]

            # Get the set of time points corresponding to DA_subgraph_i
            # This is essentially the fifteen minute points that correspond to the hour in the DA-UC problem
            ST_subgraphs_i = getsubgraphs(ST_subgraph)[(1 + (i - 1) * 4):(i * 4)]

            # Iterate over DA generators
            for k in 1:num_DA_gen
                gen_name = gen_DA[k, "Generator Name"]
                bus_name = gen_DA[k, "Bus of Connection"]
                max_ramp_up = gen_DA[k, "Max Ramp Up (MW/min)"] * 60 * dt
                max_ramp_down = gen_DA[k, "Max Ramp Down (MW/min)"] * 60 * dt

                max_startup = gen_DA[k, "Max Capacity (MW)"]
                max_shutdown = gen_DA[k, "Max Capacity (MW)"]

                gen_max_cap = gen_DA[k, "Max Capacity (MW)"]
                gen_min_cap = gen_DA[k, "Min Stable Level (MW)"]

                # Link DA generator amoutns between DA-UC and ST-UC layers
                # Ensures that the solutions don't deviate too much between layers
                for l in 1:4
                    @linkconstraint(ST_subgraph, (ST_subgraphs_i[l][:bus][bus_name][:Gen_p][gen_name] + ST_subgraphs_i[l][:bus][bus_name][:Gen_m][gen_name])
                                            - (value(DA_subgraph_i[:bus][bus_name][:Gen_p][gen_name]) + value(DA_subgraph_i[:bus][bus_name][:Gen_m][gen_name])) <= max_ramp_up / dt
                    )
                    @linkconstraint(ST_subgraph, - (ST_subgraphs_i[l][:bus][bus_name][:Gen_p][gen_name] + ST_subgraphs_i[l][:bus][bus_name][:Gen_m][gen_name])
                                            + (value(DA_subgraph_i[:bus][bus_name][:Gen_p][gen_name]) + value(DA_subgraph_i[:bus][bus_name][:Gen_m][gen_name])) <= max_ramp_up / dt
                    )
                end

                # Add ramping constraints; these constraints only apply to ST-subgraphs_i 2-4
                # For ST_subgraphs_i number 1, you have to consider the binary variables of the DA-UC level
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

                # If this is not the first time point of the day for the STUC subproblem, then link the last STUC time point
                # corresponding to the index - 1 DAUC time point to the first STUC time point corresponding to
                # the index DAUC time point. This allows for turn off and turn on conditions of the ramping constraints
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

                # Set a capacity constraint based on the generator decisions in upper layer
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

"""
    link_between_ST(graph, ST_subgraph, ST_graphs, i; dt = 0.25, day_num = 1, graph_last = nothing)

Adds linking constraints between the STUC subproblems; queries the values of previous solutions
 * graph - the full OptiGraph corresponding to a single day
 * ST_subgraph - the ST subgraph of interest to which constraints will be added
 * ST_graphs - the previous set of ST subgraphs
 * i - which ST subproblem is being passed (1-8)
 * dt - time discretization in hours
 * day_num - day number (out of 32)
 * graph_last - the full OoptiGraph for the previous day
"""
function link_between_ST_sol(graph, ST_subgraph, ST_graphs, i; dt = .25, day_num = 1, graph_last = nothing)

    """
        get_s_list(ST_graphs, bus_name, gen_name, i)

    Returns a vector of the s variables (or the value of the s variables) for several STUC subproblems; used
    for the up-time constraint; s variables are sorted in opposite direction of time point (e.g., 4:0.25:0.25)
    """
    function get_s_list(ST_graphs, bus_name, gen_name, i)
        # Builds vector of s variables (or values, in the case of previous, solved graphs)
        # depending on i, where i is the length of ST_graphs
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

    """
        get_z_list(ST_graphs, bus_name, gen_name, i)

    Returns a vector of the z variables (or the value of the z variables) for several STUC subproblems; used
    for the down-time constraint z variables are sorted in opposite direction of time point (e.g., 4:0.25:0.25)
    """
    function get_z_list(ST_graphs, bus_name, gen_name, i)
        # Builds vector of z variables (or values, in the case of previous, solved graphs)
        # depending on i, where i is the length of ST_graphs
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

    """
        link_subgraphs_uptime(graph, ST_graphs, min_up_time, bus_name, gen_name, s_list, i)

    Adds up-time linking constraints to the full, 1-day graph (graph). ST_graphs is the set of subproblem graphs (also
    used in get_s_list). Adds linking constraints for a single generator. i < 4 only applies for the first time point
    """
    function link_subgraphs_uptime(graph, ST_graphs, min_up_time, bus_name, gen_name, s_list, i)
        # The if statement is mainly used for initializing the problem. No ST-UC generator has an up time greater than 8,
        # so we only ever need 4 STUC subproblems for the up/down time linking. The first time points only link as far
        # back as there are subproblems (e.g., the third subproblem only links to the first subproblem since there are
        # no variables before this). Thus, i =1, 2, or 3 are only used for STUC subproblems 1, 2, and 3 of day 1.
        if i == 1
            # If the up time is greater than 4, then link every variable in the current subproblem (except for the first variable)
            if min_up_time >= 4
                for j in 1:15
                    subgraph_index = 17 - j
                    s_range = (j):16
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            else
                # otherwise link existing subproblems up to the current up time
                # s_list is in reverse order of time points
                # the second loop has a smaller range because it corresponds to the first variables of the subproblem
                # and so has less variables to link

                # up_time_diff is the difference between the number of time points in the subproblem and the number of time
                # points in the min_up_time
                up_time_diff = Int(4 * (4 - min_up_time))
                for j in 1:up_time_diff
                    s_range = (j):(Int(min_up_time * 4 + j - 1))
                    subgraph_index = 17 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
                for j in (up_time_diff + 1):15
                    s_range = (j):16
                    subgraph_index = 17 - j
                    # order of s_list is backwards, which is why s_range is opposite of subgraph_index
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(ST_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
            end
        elseif i == 2 # s_list is length 28
            # if min up time is between 3 and 7, the former subproblem will be needed
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
            # if the min up time is greater than 7, we need to link variables of both first and second subproblems
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

    """
        link_subgraphs_downtime(graph, ST_graphs, min_down_time, bus_name, gen_name, z_list, i)

    Adds down-time linking constraints to the full, 1-day graph (graph). ST_graphs is the set of subproblem graphs (also
    used in get_s_list). Adds linking constraints for a single generator. i < 4 only applies for the first time point
    """
    function link_subgraphs_downtime(graph, ST_graphs, min_down_time, bus_name, gen_name, z_list, i)
        # This function follows a similar format as the link_subgraphs_uptime; see the comments in that function for
        # details on the if statements
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

    # Loop through the ST generators and add ramping constraints between subproblems
    # Ramping constraints for the first time point of the subproblem uses the previous
    # time point's HAED value for its previous generator amount
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

        # Get s and z lists
        s_list = get_s_list(ST_graphs, bus_name, gen_name, ST_graphs_len)
        z_list = get_z_list(ST_graphs, bus_name, gen_name, ST_graphs_len)

        # Add up and down time constraints
        if !(isequal(0, min_up_time))
            link_subgraphs_uptime(ST_subgraph, ST_graphs, min_up_time, bus_name, gen_name, s_list, ST_graphs_len)
        end
        if !(isequal(0, min_down_time))
            link_subgraphs_downtime(ST_subgraph, ST_graphs, min_down_time, bus_name, gen_name, z_list, ST_graphs_len)
        end

        if ((day_num != 1) || (i != 1))
            # Link solutions of previous HAED layer/STUC layer with current STUC layer for the first time point of the STUC graph
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

            # Add linking constraint for binary variables using last STUC layers solution
            @linkconstraint(ST_subgraph, bus_next[:x][gen_name] - value(bus_last[:x][gen_name]) == bus_next[:s][gen_name] - bus_next[:z][gen_name])

            # Add ramping constraints
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
            # Link solutions of previous HAED layer and uppper DAUC layer with current STUC layer
            # for the first time point of the STUC graph
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

            # Define subgraphs that will be used for querying values and setting constraints
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

            # add ramping constraints
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

"""
    link_sol_over_HA(graph, tau; dt = .25)

Links the previous DAUC and STUC layer solutions to the current HA subgraph
 * graph - full day graph with 105 subgraphs
 * tau - number of HA graph to be linked (integer from 1-96)
 # dt - time discretization
"""
function link_sol_over_HA(graph, tau = 1; dt = .25)

    # Get the ST layer that the tau HA subgraph corresponds to ST_layer
    # will be 1-8 (e.g., HA subgraph 18 is on the 2nd STUC subgraph)
    ST_layer = Int(ceil(tau / 12))

    # For the STUC subgraph, this is the index that the HA_subgraph corresponds to
    # for HA subgraph (tau) = 18, we want the 2nd STUC subgraph and the start index of 6
    ST_start_index = tau%12

    # Correct the index from the % operator
    if ST_start_index == 0
        ST_start_index = 12
    end

    # get subgraphs of corresponding problems
    DA_subgraph = getsubgraphs(graph)[1]
    ST_subgraph = getsubgraphs(graph)[1 + ST_layer]
    HA_subgraph = getsubgraphs(graph)[9 + tau]

    # Loop through the STUC subproblem time points
    tau_range = ST_start_index:(ST_start_index + 4)
    for (i, index) in enumerate(tau_range)

        # Get the corresponding time point for the DAUC graph
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

            max_ramp_up = gen_DA[k, "Max Ramp Up (MW/min)"] * 60 * dt
            max_ramp_down = gen_DA[k, "Max Ramp Down (MW/min)"] * 60 * dt

            max_startup = gen_DA[k, "Max Capacity (MW)"]
            max_shutdown = gen_DA[k, "Max Capacity (MW)"]

            # Link the DAUC solutions to the HAED so that the HAED do not deviate too much from the initial decisions
            if !(isequal(0, max_ramp_up))
                @linkconstraint(HA_subgraph, (value(DA_subgraph_i[:bus][bus_name][:Gen_p][gen_name]) + value(DA_subgraph_i[:bus][bus_name][:Gen_m][gen_name])) -
                                       (HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up / dt
                )
                @linkconstraint(HA_subgraph, -(value(DA_subgraph_i[:bus][bus_name][:Gen_p][gen_name]) + value(DA_subgraph_i[:bus][bus_name][:Gen_m][gen_name])) +
                                       (HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up / dt
                )
            end

            # Set the min and max generation amounts for the generators based on the commitment decisons of the DAUC layer
            @linkconstraint(HA_subgraph, gen_max_cap * value(DA_subgraph_i[:bus][bus_name][:x][gen_name]) >=
                                   HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]
            )
            @linkconstraint(HA_subgraph, gen_min_cap * value(DA_subgraph_i[:bus][bus_name][:x][gen_name]) <=
                                   HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]
            )

            # If this is not the first time point of the full day OptiGraph, then run this code
            # If it is the first time point of the full day, then use the link_sol_over_HA_new_day instead
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

                if (tau + i - 1) % 4 == 1
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
                    # Else, the ramping constraints do not consider turn on and shutoff
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

            # Link STUC solutions for ST generators to the HA subproblem so that the solutions do not deviate by too much
            if !(isequal(0, max_ramp_up))
                @linkconstraint(HA_subgraph, (value(ST_subgraph_i[:bus][bus_name][:Gen_p][gen_name]) + value(ST_subgraph_i[:bus][bus_name][:Gen_m][gen_name])) -
                                       (HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up# / dt
                )
                @linkconstraint(HA_subgraph, -(value(ST_subgraph_i[:bus][bus_name][:Gen_p][gen_name]) + value(ST_subgraph_i[:bus][bus_name][:Gen_m][gen_name])) +
                                       (HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up# / dt
                )
            end

            # Set min and max generation limits based on commitment decisions of previous time point
            @linkconstraint(HA_subgraph, gen_max_cap * value(ST_subgraph_i[:bus][bus_name][:x][gen_name]) >=
                                   HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]
            )
            @linkconstraint(HA_subgraph, gen_min_cap * value(ST_subgraph_i[:bus][bus_name][:x][gen_name]) <=
                                   HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]
            )

            # If it is not the first time point, then set the ramping constraints
            # (If it is the first time point of the first HAED problem, the constraint is set separately
            # by the function link_sol_over_HA_new_day)
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
                if index != 1
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
            # If tau == 1, then only add ramping constraints after the first time point.
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

"""
    link_between_DA(DA_graph, DA_graph_set, day_num; dt = 1, HA_graph = nothing)

Adds up and downtime link constraints between previous the solution of previous DA problems
and the current DA problem. Also adds a ramping constraint between the current DA problem
and the final time point of the last day's DA problem and HA subproblem.
 * DA_graph - current DA graph receiving the constraints
 * DA_graph_set - set of DA subgraphs with up to the last two days graphs
 * day_num - index of day (1-32)
 * dt - time discretization for DA level
 * HA_graph - previous day's final HA graph

See also the function `link_between_ST_sol`
"""
function link_between_DA(DA_graph, DA_graph_set, day_num; dt = 1, HA_graph = nothing)

    """
        get_s_list(DA_graph_set, bus_name, gen_name, day_num)

    Returns a vector of the s variables (or the value of the s variables) for DAUC subproblems; used
    for the up-time constraint; s variables are sorted in opposite direction of time point (e.g., 4:0.25:0.25)
    """
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

    """
        get_z_list(DA_graph_set, bus_name, gen_name, i)

    Returns a vector of the z variables (or the value of the z variables) for DAUC subproblems; used
    for the down-time constraint z variables are sorted in opposite direction of time point (e.g., 4:0.25:0.25)
    """
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

    """
        link_subgraphs_uptime(graph, DA_graphs, min_up_time, bus_name, gen_name, s_list, i)

    Adds up-time linking constraints to the DA subproblem graphs. DA_graphs is the set of subproblem graphs (also
    used in get_s_list). Adds linking constraints for a single generator.
    """
    function link_subgraphs_uptime(graph, DA_graphs, min_up_time, bus_name, gen_name, s_list, i)
        # If i = 1, it is the first day, and we must account for their not being a previous graph to pull solutions from
        if i == 1
            if min_up_time > 24
                for j in 1:24
                    subgraph_index = 26 - j
                    s_range = (j):25
                    @linkconstraint(graph, sum(s_list[s_range]) <= getsubgraphs(DA_graphs[1])[subgraph_index][:bus][bus_name][:x][gen_name])
                end
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
        # if i == 2, it is the second day. If the min_up time is greater than 24, we must account for there
        # being fewer time points in the problem than the length of the minimum up or down time
        elseif i == 2
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

    """
        link_subgraphs_downtime(graph, DA_graphs, min_down_time, bus_name, gen_name, z_list, i)

    Adds down-time linking constraints to the DA subproblem graphs. DA_graphs is the set of subproblem graphs (also
    used in get_s_list). Adds linking constraints for a single generator.
    """
    function link_subgraphs_downtime(graph, DA_graphs, min_down_time, bus_name, gen_name, z_list, i)
        # If i = 1, it is the first day, and we must account for their not being a previous graph to pull solutions from
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
        # if i == 2, it is the second day. If the min_up time is greater than 24, we must account for there
        # being fewer time points in the problem than the length of the minimum up or down time
        elseif i == 2
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
    # Add constraints to each DA generator
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

            # Add constraints on binary variables
            @linkconstraint(DA_graph, bus_next[:x][gen_name] - value(bus_last[:x][gen_name]) == bus_next[:s][gen_name] - bus_next[:z][gen_name])

            # Add ramping constraints
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

"""
    link_sol_over_HA_new_day(graph, graph_last, dt = .25, monolithic = true)

Links the previous day's solutions to the 1st HA subgraph of the day
"""
function link_sol_over_HA_new_day(graph, graph_last; dt = .25, monolithic = true)

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

        max_ramp_up = gen_DA[k, "Max Ramp Up (MW/min)"] * 60 * dt
        max_ramp_down = gen_DA[k, "Max Ramp Down (MW/min)"] * 60 * dt

        max_startup = gen_DA[k, "Max Capacity (MW)"]
        max_shutdown = gen_DA[k, "Max Capacity (MW)"]

        bus_next = getsubgraphs(HA_subgraph_next)[1][:bus][bus_name]
        bus_last = getsubgraphs(HA_subgraph_last)[1][:bus][bus_name]
        bus_last_gen_p = value(bus_last[:Gen_p][gen_name])
        bus_last_gen_m = value(bus_last[:Gen_m][gen_name])

        bus_next_DA = getsubgraphs(DA_subgraph_next)[1][:bus][bus_name]
        bus_last_DA = getsubgraphs(DA_subgraph_last)[24][:bus][bus_name]

        # If this is the monolithic problem, then some DA variables
        # will not have a solution on them yet
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

        # Add ramping constraints using the previous day's DA decisions,
        # the current day's DA decision, and the previous day's HA level
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

        # If this is the monolithic problem, then some ST variables
        # will not have a solution on them yet
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

        # Add ramping constraints using the previous day's ST decisions,
        # the current day's ST decision, and the previous day's HA level
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
