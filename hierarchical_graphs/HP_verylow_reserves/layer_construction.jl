"""
    build_bus_nodes(graph, tau, D, xi, dt, gen_com_set, gen, gen_conv_set; HA::Bool = HA)

Builds the bus system for a single time point of any of the DAUC, STUC, or HA subproblems.
Takes the following arguments:
 * graph - empty OptiGraph for the given time point
 * tau - time point of the subproblem (used for indexing the renewables)
 * D - demand at each bus (has length 118)
 * xi - dictionary of the renewable generation; keys are renewable generation names
 * dt - discretization size (in hours)
 * gen_com_set - set of generators to be committed in subproblem (for HA layer, this must be empty)
 * gen - set of generators in the subproblem (DA-UC only contains DA-UC and renewable generators)
 * gen_conv_set - set of conventional generators in the subproblem
 * HA - indicator whether this is the HA-ED layer or not; if it is, the objective is different than the other layers
"""
function build_bus_nodes(graph, tau, D, xi, dt, gen_com_set, gen, gen_conv_set; HA::Bool = HA)
    # Added nodes for each bus and for each line
    @optinode(graph, bus[bus_data[:, "Bus Name"]])
    @optinode(graph, line[line_data[:, "Line Name"]])

    for (i, node) in enumerate(bus)
        # Get the set of generators at the given bus
        gen_bus_bool = gen[:, "Bus of Connection"] .== bus_data[i, "Bus Name"]

        gen_bus = gen[gen_bus_bool, "Generator Name"]

        # Define power to grid (Gen_p) and curtailed/overgenerated power (Gen_m)
        @variable(node, 0 <= Gen_p[gen_bus])
        @variable(node, 0 <= Gen_m[gen_bus])

        # Set upper bounds; this is added to help the MIP solver by bounding the variables
        # Generator maximum capacity is a constraint set later
        for gen_name in gen_bus
            gen_index = findfirst(x -> x == gen_name, gen_data[:, "Generator Name"])
            max_cap = gen_data[gen_index, "Max Capacity (MW)"]
            set_upper_bound(Gen_p[gen_name], max_cap)
            set_upper_bound(Gen_m[gen_name], max_cap)
        end

        # Get the set of generators at this bus which are to be committed
        gen_com_bool = findall(x -> (x in gen_com_set[:, "Generator Name"]), gen_bus)
        gen_com = gen_bus[gen_com_bool]

        # Define binary variables
        # For HA-ED layer, gen_com is empty, so these variables won't exist
        @variable(node, x[gen_com], Bin)
        @variable(node, s[gen_com], Bin)
        @variable(node, z[gen_com], Bin)

        gen_com_max_cap = gen[gen_bus_bool, "Max Capacity (MW)"][gen_com_bool]
        gen_com_min_cap = gen[gen_bus_bool, "Min Stable Level (MW)"][gen_com_bool]
        gen_com_start_cost = gen[gen_bus_bool, "Start Cost (\$)"][gen_com_bool]
        gen_com_no_load_cost = gen[gen_bus_bool, "No-Load Cost (\$/hr)"][gen_com_bool]

        # Define the maximum/minimum capacity of the generator; dependent on whether the generator is on or not
        for j in 1:length(gen_com)
            @constraint(node, gen_com_max_cap[j] * x[gen_com[j]] >= Gen_p[gen_com[j]] + Gen_m[gen_com[j]])
            @constraint(node, gen_com_min_cap[j] * x[gen_com[j]] <= Gen_p[gen_com[j]] + Gen_m[gen_com[j]])
        end

        # Get the renewable generators at the given bus
        gen_ren_bus_bool = findall(x -> (x in gen_r[:, "Generator Name"]), gen_bus)
        gen_ren_bus = gen_bus[gen_ren_bus_bool]

        # Fix the generation of the renewable generators
        for j in 1:length(gen_ren_bus)
            @constraint(node, Gen_p[gen_ren_bus[j]] + Gen_m[gen_ren_bus[j]] == xi[gen_ren_bus[j]][tau])
        end

        # Define bus variables
        @variable(node, -pi <= theta <= pi)
        @variable(node, 0 <= D_shed <= 20000)

        gen_conv_bus_bool = findall(x -> x in gen_conv_set[:, "Generator Name"], gen_bus)
        gen_conv_bus = gen_bus[gen_conv_bus_bool]

        # Define objective function
        if HA
            var_gen_cost = gen[gen_bus_bool, "Variable Cost (\$/MWh)"]

            @objective(node, Min,
                sum(var_gen_cost[k] * Gen_p[gen_bus[k]] * dt for k in 1:length(gen_bus)) +
                sum(phi_o * Gen_m[gen_conv_bus[k]] for k in 1:length(gen_conv_bus)) * dt +
                sum(phi_c * Gen_m[gen_ren_bus[k]] for k in 1:length(gen_ren_bus)) * dt +
                phi_u * D_shed * dt
            )
        else
            var_gen_cost = gen[gen_bus_bool, "Variable Cost (\$/MWh)"][gen_com_bool]

            @objective(node, Min,
                sum(gen_com_start_cost[k] * s[gen_com[k]] + gen_com_no_load_cost[k] * x[gen_com[k]] * dt +
                var_gen_cost[k] * Gen_p[gen_com[k]] * dt for k in 1:length(gen_com)) +
                sum(phi_o * Gen_m[gen_com[k]] for k in 1:length(gen_com)) * dt +
                sum(phi_c * Gen_m[gen_ren_bus[k]] for k in 1:length(gen_ren_bus)) * dt +
                phi_u * D_shed * dt
            )
        end
    end

    # Add constraints for minimum and maximum flow on each transmission line
    for (j, node) in enumerate(line)
        Fmin = line_data[j, "Min Flow (MW)"]
        Fmax = line_data[j, "Max Flow (MW)"]
        @variable(node, Fmin <= F <= Fmax)

        from_bus = line_data[j, "Bus from"]
        to_bus = line_data[j, "Bus to"]
        @linkconstraint(graph, F == B[j] * (bus[from_bus][:theta] - bus[to_bus][:theta]))
    end

    # Constrain flow at each bus
    for j in 1:length(bus)
        bus_name = bus_data[j, "Bus Name"]

        bus_in_set = findall(x -> x == bus_name, line_data[:, "Bus to"]) # set of all bus nodes that go to current bus
        bus_out_set = findall(x -> x == bus_name, line_data[:, "Bus from"]) # set of all bus nodes that go from current node
        @linkconstraint(graph,
            sum(line[line_names][:F] for line_names in line_data[bus_in_set, "Line Name"]) -
            sum(line[line_names][:F] for line_names in line_data[bus_out_set, "Line Name"]) +
            sum(k for k in bus[bus_name][:Gen_p]) +
            bus[bus_name][:D_shed]
            == D[j]
        )
    end
end

"""
    build_bus_over_time(graph, t_len, D, xi, dt, gen_com_set, gen, gen_conv_set; offset = 0, HA::Bool = false)

Builds a subgraph for each time point of a given subproblem (e.g., for the DA-UC problem, it builds 25
time point graphs and adds them to the subproblem graph). Takes the following arguments:
 * graph - OptiGraph of the subproblem
 * t_len - the number of time points to build a graph for
 * D - DataFrame of demands for each bus for each time point; if t_len is 25, this should be a DataFrame size 118 x 25
 * xi - the set of renewable generation data for the subproblem
 * dt - discretization size in hours
 * gen_com_set - set of generators to be committed in subproblem (for HA layer, this must be empty)
 * gen - set of generators in the subproblem (DA-UC only contains DA-UC and renewable generators)
 * gen_conv_set - set of conventional generators in the subproblem
 * offset - offset used for getting the correct index from xi; for example, if this is the second day of the DA-UC
   problem, then this should be set to 24
 * HA - indicator whether this is the HA-ED layer or not; if it is, the objective is different than the other layers
"""
function build_bus_over_time(graph, t_len, D, xi, dt, gen_com_set, gen, gen_conv_set; offset = 0, HA::Bool = false)
    # Iterate over t_len
    for i in 1:t_len
        # Define OptiGraph for the time point
        t_graph = OptiGraph()
        # Add the necessary data to the OptiGraph
        build_bus_nodes(t_graph, i + offset, D[:, i], xi, dt, gen_com_set, gen, gen_conv_set, HA = HA)
        # Add the time point graph to the subproblem graph
        add_subgraph!(graph, t_graph)
    end
end

"""
    link_DA_over_time(graph, gen)

Add all necessary linking constraints within the DA-UC subproblem graph
 * graph - The subproblem graph for the DA-UC subproblem
 * gen - The set of all generators that show up in the DA-UC subproblem (should be DA-UC committed generators + renewables)
"""
function link_DA_over_time(graph, gen)
    # Get all subgraphs (each a time point) in the DAUC problem
    DAUC_subgraphs = all_subgraphs(graph)

    # Get the number of subgraphs (number of time points)
    t_len = length(DAUC_subgraphs)

    # Iterate through t_len - 1
    for i in 1:(t_len - 1)
        # Iterate through the number of buses (118)
        for j in 1:length(DAUC_subgraphs[1][:bus])
            # Get the bus name and the set of generators on that bus
            bus_name = bus_data[j, "Bus Name"]
            gen_set_bool = gen[:, "Bus of Connection"] .== bus_name

            gen_set = gen[gen_set_bool, "Generator Name"]
            max_ramp_up = gen[gen_set_bool, "Max Ramp Up (MW/min)"]
            max_ramp_down = gen[gen_set_bool, "Max Ramp Down (MW/min)"]
            max_cap = gen[gen_set_bool, "Max Capacity (MW)"]
            min_cap = gen[gen_set_bool, "Min Stable Level (MW)"]

            max_startup = gen[gen_set_bool, "Max Capacity (MW)"]
            max_shutdown = gen[gen_set_bool, "Max Capacity (MW)"]

            # Get the set of conventional generators (i.e., remove renewables)
            gen_conv_bool = findall(x -> (x in gen_DA[:, "Generator Name"]) || (x in gen_ST[:, "Generator Name"]), gen_set)
            gen_conv = gen_set[gen_conv_bool]
            max_ramp_up_conv = max_ramp_up[gen_conv_bool] .* 60 #per hour; max_ramp_up is in minutes
            max_ramp_down_conv = max_ramp_down[gen_conv_bool] .* 60 #per hour
            max_cap_conv = max_cap[gen_conv_bool]
            min_cap_conv = min_cap[gen_conv_bool]

            max_startup_conv = max_startup[gen_conv_bool]
            max_shutdown_conv = max_startup[gen_conv_bool]

            bus_next = DAUC_subgraphs[i + 1][:bus][bus_name]
            bus_last = DAUC_subgraphs[i][:bus][bus_name]

            # Add linking constraints between the time point graphs
            for k in 1:length(gen_conv)
                # Constrain x, s, and z variables
                @linkconstraint(graph, bus_next[:x][gen_conv[k]] - bus_last[:x][gen_conv[k]] ==
                    bus_next[:s][gen_conv[k]] - bus_next[:z][gen_conv[k]]
                )

                # Define ramping constraints
                @linkconstraint(graph, bus_next[:Gen_p][gen_conv[k]] + bus_next[:Gen_m][gen_conv[k]] -
                    bus_last[:Gen_p][gen_conv[k]] - bus_last[:Gen_m][gen_conv[k]] <=
                    (max_startup_conv[k] - max_ramp_up_conv[k] - min_cap_conv[k]) * bus_next[:s][gen_conv[k]] +
                    (max_ramp_up_conv[k] + min_cap_conv[k]) * bus_next[:x][gen_conv[k]] - min_cap_conv[k] * bus_last[:x][gen_conv[k]]
                )

                @linkconstraint(graph, bus_last[:Gen_p][gen_conv[k]] + bus_last[:Gen_m][gen_conv[k]] -
                    bus_next[:Gen_p][gen_conv[k]] - bus_next[:Gen_m][gen_conv[k]] <=
                    (max_shutdown_conv[k] - max_ramp_down_conv[k] - min_cap_conv[k]) * bus_next[:z][gen_conv[k]] +
                    (max_ramp_down_conv[k] + min_cap_conv[k]) * bus_last[:x][gen_conv[k]] - min_cap_conv[k] * bus_next[:x][gen_conv[k]]
                )
            end
        end
    end
end

"""
    link_ST_over_time(graph, dt = .25)

Add all necessary linking constraints within the ST-UC subproblem graph
 * graph - The subproblem graph for the ST-UC subproblem
 * dt - discretization size in hours
"""
function link_ST_over_time(graph, dt = .25)

    STUC_subgraphs = all_subgraphs(graph)

    t_len = length(STUC_subgraphs)

    for i in 1:(t_len - 1)
        for j in 1:length(STUC_subgraphs[1][:bus])
            bus_name = bus_data[j, "Bus Name"]
            gen_set_bool = gen_data[:, "Bus of Connection"] .== bus_name

            gen_set = gen_data[gen_set_bool, "Generator Name"]
            max_ramp_up = gen_data[gen_set_bool, "Max Ramp Up (MW/min)"]
            max_ramp_down = gen_data[gen_set_bool, "Max Ramp Down (MW/min)"]
            max_cap = gen_data[gen_set_bool, "Max Capacity (MW)"]
            min_cap = gen_data[gen_set_bool, "Min Stable Level (MW)"]

            max_startup = gen_data[gen_set_bool, "Max Capacity (MW)"]
            max_shutdown = gen_data[gen_set_bool, "Max Capacity (MW)"]

            # conventional set
            gen_ST_bool = findall(x -> (x in gen_ST[:, "Generator Name"]), gen_set)
            gen_conv = gen_set[gen_ST_bool]
            max_ramp_up_conv = max_ramp_up[gen_ST_bool] .* 60 .* dt
            max_ramp_down_conv = max_ramp_down[gen_ST_bool] .* 60 .* dt
            max_cap_conv = max_cap[gen_ST_bool]
            min_cap_conv = min_cap[gen_ST_bool]

            max_startup_conv = max_startup[gen_ST_bool]
            max_shutdown_conv = max_shutdown[gen_ST_bool]

            bus_next = STUC_subgraphs[i + 1][:bus][bus_name]
            bus_last = STUC_subgraphs[i][:bus][bus_name]

            for k in 1:length(gen_conv)
                if gen_conv[k] in gen_ST[:, "Generator Name"]
                    @linkconstraint(graph, bus_next[:x][gen_conv[k]] - bus_last[:x][gen_conv[k]] ==
                        bus_next[:s][gen_conv[k]] - bus_next[:z][gen_conv[k]]
                    )
                end

                @linkconstraint(graph, bus_next[:Gen_p][gen_conv[k]] + bus_next[:Gen_m][gen_conv[k]] -
                    bus_last[:Gen_p][gen_conv[k]] - bus_last[:Gen_m][gen_conv[k]] <=
                    (max_startup_conv[k] - max_ramp_up_conv[k] - min_cap_conv[k]) * bus_next[:s][gen_conv[k]] +
                    (max_ramp_up_conv[k] + min_cap_conv[k]) * bus_next[:x][gen_conv[k]] - min_cap_conv[k] * bus_last[:x][gen_conv[k]]
                )

                @linkconstraint(graph, bus_last[:Gen_p][gen_conv[k]] + bus_last[:Gen_m][gen_conv[k]] -
                    bus_next[:Gen_p][gen_conv[k]] - bus_next[:Gen_m][gen_conv[k]] <=
                    (max_shutdown_conv[k] - max_ramp_down_conv[k] - min_cap_conv[k]) * bus_next[:z][gen_conv[k]] +
                    (max_ramp_down_conv[k] + min_cap_conv[k]) * bus_last[:x][gen_conv[k]] - min_cap_conv[k] * bus_next[:x][gen_conv[k]]
                )

            end
        end
    end

    num_ST_gen = size(gen_ST, 1)
end

"""
    link_DA_to_ST(graph::OptiGraph, tau_0 = 1:4, ST_sub_num = 1, dt = .25, gen_com_set = gen_ST)

Links the variables in the DA level to the ST level.
 * graph - full day OptiGraph with all three layers (should have DAUC as first OptiGraph, STUC as next 8)
 * tau_0 - range of hours in the DA level that ST-UC solves for (e.g., 1-4, 4-7, etc.)
 * ST_sub_num - number of ST_graph (e.g., for 1-4, should be 1, for 4-7, should be 2, etc.)
 * dt - discretization size (in hours)
 * gen_com_set - the set of generators to be committed
"""
function link_DA_to_ST(graph::OptiGraph, tau_0 = 1:4, ST_sub_num = 1, dt = .25, gen_com_set = gen_ST)

    # Get the DAUC subgraph
    DA_subgraph = getsubgraphs(graph)[1]
    # Get the STUC subgraph that needs to be linked
    ST_subgraph = getsubgraphs(graph)[1 + ST_sub_num]

    # get the number of DA generators
    num_DA_gen = size(gen_DA, 1)

    # Iterate through the DA-UC time points that correspond to the problem
    for (i, index) in enumerate(tau_0)

        if index <= length(getsubgraphs(DA_subgraph))
            # Get the time point subgraph from DA_subgraph
            DA_subgraph_i = getsubgraphs(DA_subgraph)[index]
            # Get the set of time points corresponding to the DA_subgraph_i
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

                # Link DA generator amounts between DA-UC and ST-UC layers
                # Ensures that the solutions don't deviate too much between layers
                for l in 1:4
                    @linkconstraint(graph, (ST_subgraphs_i[l][:bus][bus_name][:Gen_p][gen_name] + ST_subgraphs_i[l][:bus][bus_name][:Gen_m][gen_name])
                                            - (DA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + DA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up / dt
                    )
                    @linkconstraint(graph, - (ST_subgraphs_i[l][:bus][bus_name][:Gen_p][gen_name] + ST_subgraphs_i[l][:bus][bus_name][:Gen_m][gen_name])
                                            + (DA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + DA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up / dt
                    )
                end

                # Add ramping constraints; these constraints only apply to ST_subgraphs_i 2-4
                # For ST_subgraphs_i number 1, you have to consider the binary variables of the DA-UC level
                for l in 2:4
                    bus_next = ST_subgraphs_i[l][:bus][bus_name]
                    bus_last = ST_subgraphs_i[l - 1][:bus][bus_name]
                    @linkconstraint(graph, (bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name]) -
                                           (bus_last[:Gen_p][gen_name] + bus_last[:Gen_m][gen_name]) <=
                                           max_ramp_up
                    )
                    @linkconstraint(graph, (bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name]) -
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
                    @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                        bus_last[:Gen_p][gen_name] - bus_last[:Gen_m][gen_name] <=
                        (max_startup - max_ramp_up - gen_min_cap) * bus_next_DA[:s][gen_name] +
                        (max_ramp_up + gen_min_cap) * bus_next_DA[:x][gen_name] - gen_min_cap * bus_last_DA[:x][gen_name]
                    )

                    @linkconstraint(graph, bus_last[:Gen_p][gen_name] + bus_last[:Gen_m][gen_name] -
                        bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                        (max_shutdown - max_ramp_down - gen_min_cap) * bus_next_DA[:z][gen_name] +
                        (max_ramp_down + gen_min_cap) * bus_last_DA[:x][gen_name] - gen_min_cap * bus_next_DA[:x][gen_name]
                    )
                end

                # Set a capacity constraint based on the generator decisions in upper layer
                for (j, subgraph) in enumerate(ST_subgraphs_i)
                    @linkconstraint(graph, gen_max_cap * DA_subgraph_i[:bus][bus_name][:x][gen_name] >=
                                           subgraph[:bus][bus_name][:Gen_p][gen_name] + subgraph[:bus][bus_name][:Gen_m][gen_name]
                    )
                    @linkconstraint(graph, gen_min_cap * DA_subgraph_i[:bus][bus_name][:x][gen_name] <=
                                           subgraph[:bus][bus_name][:Gen_p][gen_name] + subgraph[:bus][bus_name][:Gen_m][gen_name]
                    )
                end
            end
        end
    end
end

"""
    link_between_ST(graph; dt = .25, day_num = 1, graph_last = nothing)

Adds linking constraints between the STUC subproblems
 * graph - the full OptiGraph corresopnding to a single day
 * dt - time discretization in hours
 * day_num - day number (out of 32); really just used to ensure the linking behaves differently on the first day (could/should be a Boolean)
 * graph_last - the full OptiGraph for the previous day
"""
function link_between_ST(graph; dt = .25, day_num = 1, graph_last = nothing)

    """
        get_s_list(ST_graphs, bus_name, gen_name, i)

    Returns a vector of the s variables for several STUC subproblems; used for the up-time constraint
    s variables are sorted in opposite direction of time point (e.g., 4:0.25:0.25, not 0.25:0.25:4)
    """
    function get_s_list(ST_graphs, bus_name, gen_name, i)
        # Builds vector of s variables depending on i, where i is the length of ST_graphs
        if i == 1
            s_list = Vector{Any}([getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1])
        elseif i == 2
            s_list1 = [getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1]
            s_list2 = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list = Vector{Any}(vcat(s_list1, s_list2))
        elseif i == 3
            s_list1 = [getsubgraphs(ST_graphs[3])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1]
            s_list2 = [getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list3 = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list = Vector{Any}(vcat(s_list1, s_list2, s_list3))
        else
            s_list1 = [getsubgraphs(ST_graphs[i])[j][:bus][bus_name][:s][gen_name] for j in 16:-1:1]
            s_list2 = [getsubgraphs(ST_graphs[i - 1])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list3 = [getsubgraphs(ST_graphs[i - 2])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list4 = [getsubgraphs(ST_graphs[i - 3])[j][:bus][bus_name][:s][gen_name] for j in 12:-1:1]
            s_list = Vector{Any}(vcat(s_list1, s_list2, s_list3, s_list4))
        end
        return s_list
    end

    """
        get_z_list(ST_graphs, bus_name, gen_name, i)

    Returns a vector of the z variables for several STUC subproblems; used for the down-time constraint
    z variables are sorted in opposite direction of time point (e.g., 4:0.25:0.25, not 0.25:0.25:4)
    """
    function get_z_list(ST_graphs, bus_name, gen_name, i)
        # builds vector of z variables depending on i, where i is the length of ST_graphs
        if i == 1
            z_list = Vector{Any}([getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1])
        elseif i == 2
            z_list1 = [getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1]
            z_list2 = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list = Vector{Any}(vcat(z_list1, z_list2))
        elseif i == 3
            z_list1 = [getsubgraphs(ST_graphs[3])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1]
            z_list2 = [getsubgraphs(ST_graphs[2])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list3 = [getsubgraphs(ST_graphs[1])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list = Vector{Any}(vcat(z_list1, z_list2, z_list3))
        else
            z_list1 = [getsubgraphs(ST_graphs[i])[j][:bus][bus_name][:z][gen_name] for j in 16:-1:1]
            z_list2 = [getsubgraphs(ST_graphs[i - 1])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list3 = [getsubgraphs(ST_graphs[i - 2])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list4 = [getsubgraphs(ST_graphs[i - 3])[j][:bus][bus_name][:z][gen_name] for j in 12:-1:1]
            z_list = Vector{Any}(vcat(z_list1, z_list2, z_list3, z_list4))
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

    num_ST_gen = size(gen_ST, 1)

    # Loop through all STUC subproblems
    for i in 1:8
        if day_num != 1
            # if the day_num is not 1, then we need to get teh values of the previous days s_variables
            if i == 1
                # set the subgraph set to include the previous days STUC subproblem graphs
                ST_subgraph_set1 = getsubgraphs(graph_last)[7:9]
                ST_subgraph_set2 = getsubgraphs(graph)[2:2]
                ST_graphs = vcat(ST_subgraph_set1, ST_subgraph_set2)
            elseif i == 2
                # set the subgraph set to include the previous days STUC subproblem graphs
                ST_subgraph_set1 = getsubgraphs(graph_last)[8:9]
                ST_subgraph_set2 = getsubgraphs(graph)[2:3]
                ST_graphs = vcat(ST_subgraph_set1, ST_subgraph_set2)
            elseif i == 3
                # set the subgraph set to include the previous days STUC subproblem graphs
                ST_subgraph_set1 = getsubgraphs(graph_last)[9:9]
                ST_subgraph_set2 = getsubgraphs(graph)[2:4]
                ST_graphs = vcat(ST_subgraph_set1, ST_subgraph_set2)
            else
                # ST_graphs only includes the current day's graphs
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

            if day_num != 1
                s_list = get_s_list(ST_graphs, bus_name, gen_name, ST_graphs_len)
                z_list = get_z_list(ST_graphs, bus_name, gen_name, ST_graphs_len)

                # We need to convert the previous day's variables to values
                # Set the index where the previous day's s and z variables start
                if i == 1
                    start_index = 17
                elseif i == 2
                    start_index = 29
                elseif i == 3
                    start_index = 41
                else
                    start_index = 53
                end
                for k in start_index:52
                    s_list[k] = value(s_list[k])
                    z_list[k] = value(z_list[k])
                end

                # If up and down times are not zero, then add constraints for up and down time
                if !(isequal(0, min_up_time))
                    link_subgraphs_uptime(graph, ST_graphs, min_up_time, bus_name, gen_name, s_list, ST_graphs_len)
                end
                if !(isequal(0, min_down_time))
                    link_subgraphs_downtime(graph, ST_graphs, min_down_time, bus_name, gen_name, z_list, ST_graphs_len)
                end
            else
                # If it's the first day, the functions are easy to call:
                s_list = get_s_list(ST_graphs, bus_name, gen_name, i)
                z_list = get_z_list(ST_graphs, bus_name, gen_name, i)
                if !(isequal(0, min_up_time))
                    link_subgraphs_uptime(graph, ST_graphs, min_up_time, bus_name, gen_name, s_list, i)
                end
                if !(isequal(0, min_down_time))
                    link_subgraphs_downtime(graph, ST_graphs, min_down_time, bus_name, gen_name, z_list, i)
                end
            end
        end

        # Loop through the ST generators and add ramping constraints between subproblems
        # Ramping constraints for the first time point of the subproblem uses the previous days
        # HAED value for its previous generator amount
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

            ST_graphs = getsubgraphs(graph)[2:9]

            # Add ramping constraints between the first time point of the ST-UC layer
            # and the previous time point's HAED layer
            if i != 1
                # if this is not the first time point, link to the previous time point's graphs (within the object graph)
                bus_next = getsubgraphs(ST_graphs[i])[1][:bus][bus_name]
                bus_last = getsubgraphs(ST_graphs[i - 1])[12][:bus][bus_name]
                HA_subgraph = getsubgraphs(graph)[9 + (i - 1) * 12]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]

                # link binary variables
                @linkconstraint(graph, bus_next[:x][gen_name] - bus_last[:x][gen_name] == bus_next[:s][gen_name] - bus_next[:z][gen_name])

                # add ramping constraints
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
                # if i == 1 and the day_num is not 1, then link to the previous graph
                bus_next = getsubgraphs(ST_graphs[1])[1][:bus][bus_name]
                bus_last = getsubgraphs(getsubgraphs(graph_last)[9])[12][:bus][bus_name]
                HA_subgraph = getsubgraphs(graph_last)[105]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]

                # link binary variables
                @linkconstraint(graph, bus_next[:x][gen_name] - value(bus_last[:x][gen_name]) == bus_next[:s][gen_name] - bus_next[:z][gen_name])

                # add ramping constraints
                @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                    value(bus_last_HA[:Gen_p][gen_name]) - value(bus_last_HA[:Gen_m][gen_name]) <=
                    (max_startup - max_ramp_up - min_cap) * bus_next[:s][gen_name] +
                    (max_ramp_up + min_cap) * bus_next[:x][gen_name] - min_cap * value(bus_last[:x][gen_name])
                )

                @linkconstraint(graph, value(bus_last_HA[:Gen_p][gen_name]) + value(bus_last_HA[:Gen_m][gen_name]) -
                    bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                    (max_shutdown - max_ramp_down - min_cap) * bus_next[:z][gen_name] +
                    (max_ramp_down + min_cap) * value(bus_last[:x][gen_name]) - min_cap * bus_next[:x][gen_name]
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

    # Loop through the DA generators and add ramping constraints between subproblems
    # Ramping constraints for the first time point of the subproblem uses the previous days
    # HAED value for its previous generator amount
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
                # Add ramping constraints between the first time point of the STUC layer
                # and the previous day's HAED layer and previous day's DAUC layer
                bus_next = getsubgraphs(ST_graphs[i])[1][:bus][bus_name]
                bus_last = getsubgraphs(getsubgraphs(graph_last)[9])[12][:bus][bus_name]

                HA_subgraph = getsubgraphs(graph_last)[105]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]

                DA_subgraph_last = getsubgraphs(graph_last)[1]
                bus_last_DA = getsubgraphs(DA_subgraph_last)[24][:bus][bus_name]

                DA_subgraph_next = getsubgraphs(graph)[1]
                bus_next_DA = getsubgraphs(DA_subgraph_next)[1][:bus][bus_name]

                # Add ramping constraints
                @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                    value(bus_last_HA[:Gen_p][gen_name]) - value(bus_last_HA[:Gen_m][gen_name]) <=
                    (max_startup - max_ramp_up / dt - min_cap) * bus_next_DA[:s][gen_name] +
                    (max_ramp_up / dt + min_cap) * bus_next_DA[:x][gen_name] - min_cap * value(bus_last_DA[:x][gen_name])
                )

                @linkconstraint(graph, value(bus_last_HA[:Gen_p][gen_name]) + value(bus_last_HA[:Gen_m][gen_name]) -
                    bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                    (max_shutdown - max_ramp_down / dt - min_cap) * bus_next_DA[:z][gen_name] +
                    (max_ramp_down / dt + min_cap) * value(bus_last_DA[:x][gen_name]) - min_cap * bus_next_DA[:x][gen_name]
                )
            else
                # Add ramping constraints and link between DAUC, STUC, and HAED layers
                bus_next = getsubgraphs(ST_graphs[i])[1][:bus][bus_name]
                bus_last = getsubgraphs(ST_graphs[i - 1])[12][:bus][bus_name]

                HA_subgraph = getsubgraphs(graph)[9 + (i - 1) * 12]
                bus_last_HA = getsubgraphs(HA_subgraph)[1][:bus][bus_name]

                DA_subgraph = getsubgraphs(graph)[1]
                bus_last_DA = getsubgraphs(DA_subgraph)[(i - 1) * 3][:bus][bus_name]

                DA_subgraph = getsubgraphs(graph)[1]
                bus_next_DA = getsubgraphs(DA_subgraph)[1 + (i - 1) * 3][:bus][bus_name]

                # add ramping constraints
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
end

"""
    link_over_HA(graph, tau = 1; dt = .25)

Adds linking constraints to the HA subgraphs that links between HAED subproblems. Also links to the
generator commitment decisions for ramping constraints. graph is the OptiGraph for the entire day.
tau is the HAED subproblem number (out of 96).
"""
function link_over_HA(graph, tau = 1; dt = .25)

    # Decides which STUC subproblem the HAED subproblem corresponds to
    ST_layer = Int(ceil(tau / 12))
    # Decides which time point of the STUC subproblem to start at
    ST_start_index = tau%12

    if ST_start_index == 0
        ST_start_index = 12
    end

    # Get subgraphs of corresponding problems
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
                @linkconstraint(graph, (DA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + DA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) -
                                       (HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up / dt
                )
                @linkconstraint(graph, -(DA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + DA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) +
                                       (HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up / dt
                )
            end

            # Set the min and max generation amounts for the generators based on the commitment decisons of the DAUC layer
            @linkconstraint(graph, gen_max_cap * DA_subgraph_i[:bus][bus_name][:x][gen_name] >=
                                   HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]
            )
            @linkconstraint(graph, gen_min_cap * DA_subgraph_i[:bus][bus_name][:x][gen_name] <=
                                   HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]
            )

            # If this is not the first time point of the full day OptiGraph, then run this code
            # If it is the first time point of the full day, then use the link_sol_over_HA_new_day instead
            if (tau != 1)
                bus_next = HA_subgraph_i[:bus][bus_name]
                if i != 1
                    bus_last = getsubgraphs(HA_subgraph)[i - 1][:bus][bus_name]
                else
                    HA_subgraph_last = getsubgraphs(graph)[9 + tau - 1]
                    bus_last = getsubgraphs(HA_subgraph_last)[1][:bus][bus_name]
                end

                # If the HAED index corresponds to the first time point of the subproblem (and not the first time point of the graph)
                # then we need to add ramping constraints that allow for turn on and shutoff.
                if tau%4 == 1
                    bus_next_DA = getsubgraphs(DA_subgraph)[DA_index][:bus][bus_name]
                    bus_last_DA = getsubgraphs(DA_subgraph)[DA_index - 1][:bus][bus_name]
                    @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                        bus_last[:Gen_p][gen_name] - bus_last[:Gen_m][gen_name] <=
                        (max_startup - max_ramp_up / dt - gen_min_cap) * bus_next_DA[:s][gen_name] +
                        (max_ramp_up / dt + gen_min_cap) * bus_next_DA[:x][gen_name] - gen_min_cap * bus_last_DA[:x][gen_name]
                    )

                    @linkconstraint(graph, bus_last[:Gen_p][gen_name] + bus_last[:Gen_m][gen_name] -
                        bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                        (max_shutdown - max_ramp_down / dt - gen_min_cap) * bus_next_DA[:z][gen_name] +
                        (max_ramp_down / dt + gen_min_cap) * bus_last_DA[:x][gen_name] - gen_min_cap * bus_next_DA[:x][gen_name]
                    )
                else
                    # Else, the ramping constraints do not consider turn on and shutoff
                    @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                        bus_last[:Gen_p][gen_name] - bus_last[:Gen_m][gen_name] <= max_ramp_up
                    )
                    @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                        bus_last[:Gen_p][gen_name] - bus_last[:Gen_m][gen_name] >= - max_ramp_down
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
                @linkconstraint(graph, (ST_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + ST_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) -
                                       (HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up# / dt
                )
                @linkconstraint(graph, -(ST_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + ST_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) +
                                       (HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]) <= max_ramp_up# / dt
                )
            end

            # Set min and max generation limits based on commitment decisions of previous time point
            @linkconstraint(graph, gen_max_cap * ST_subgraph_i[:bus][bus_name][:x][gen_name] >=
                                   HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]
            )
            @linkconstraint(graph, gen_min_cap * ST_subgraph_i[:bus][bus_name][:x][gen_name] <=
                                   HA_subgraph_i[:bus][bus_name][:Gen_p][gen_name] + HA_subgraph_i[:bus][bus_name][:Gen_m][gen_name]
            )

            # If it is not the first time point, then set the ramping constraints
            # (If it is the first time point of the first HAED problem, the constraint is set separately
            # by the function link_sol_over_HA_new_day)
            if (tau != 1)
                bus_next = HA_subgraph_i[:bus][bus_name]
                if i != 1
                    bus_last = getsubgraphs(HA_subgraph)[i - 1][:bus][bus_name]
                else
                    HA_subgraph_last = getsubgraphs(graph)[9 + tau - 1]
                    bus_last = getsubgraphs(HA_subgraph_last)[1][:bus][bus_name]
                end
                bus_next_ST = getsubgraphs(ST_subgraph)[index][:bus][bus_name]
                if ST_start_index != 1
                    bus_last_ST = getsubgraphs(ST_subgraph)[index - 1][:bus][bus_name]
                else
                    ST_subgraph_last = getsubgraphs(graph)[1 + ST_layer - 1]
                    bus_last_ST = getsubgraphs(ST_subgraph_last)[12][:bus][bus_name]
                end
                @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                    bus_last[:Gen_p][gen_name] - bus_last[:Gen_m][gen_name] <=
                    (max_startup - max_ramp_up - gen_min_cap) * bus_next_ST[:s][gen_name] +
                    (max_ramp_up + gen_min_cap) * bus_next_ST[:x][gen_name] - gen_min_cap * bus_last_ST[:x][gen_name]
                )

                @linkconstraint(graph, bus_last[:Gen_p][gen_name] + bus_last[:Gen_m][gen_name] -
                    bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                    (max_shutdown - max_ramp_down - gen_min_cap) * bus_next_ST[:z][gen_name] +
                    (max_ramp_down + gen_min_cap) * bus_last_ST[:x][gen_name] - gen_min_cap * bus_next_ST[:x][gen_name]
                )
            # If tau == 1, then only add ramping constraints after the first time point.
            elseif i != 1
                bus_next = HA_subgraph_i[:bus][bus_name]
                bus_last = getsubgraphs(HA_subgraph)[i - 1][:bus][bus_name]
                bus_last_gen_p = bus_last[:Gen_p][gen_name]
                bus_last_gen_m = bus_last[:Gen_m][gen_name]

                bus_next_ST = getsubgraphs(ST_subgraph)[index][:bus][bus_name]
                bus_last_ST = getsubgraphs(ST_subgraph)[index - 1][:bus][bus_name]
                @linkconstraint(graph, bus_next[:Gen_p][gen_name] + bus_next[:Gen_m][gen_name] -
                    bus_last_gen_p - bus_last_gen_m <=
                    (max_startup - max_ramp_up - gen_min_cap) * bus_next_ST[:s][gen_name] +
                    (max_ramp_up + gen_min_cap) * bus_next_ST[:x][gen_name] - gen_min_cap * bus_last_ST[:x][gen_name]
                )

                @linkconstraint(graph, bus_last_gen_p + bus_last_gen_m -
                    bus_next[:Gen_p][gen_name] - bus_next[:Gen_m][gen_name] <=
                    (max_shutdown - max_ramp_down - gen_min_cap) * bus_next_ST[:z][gen_name] +
                    (max_ramp_down + gen_min_cap) * bus_last_ST[:x][gen_name] - gen_min_cap * bus_next_ST[:x][gen_name]
                )
            end
        end
    end
end
