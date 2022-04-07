using JLD2
using Plasmo
# Define a function that will create the pipeline for a single scenario
function make_stoch_subgraph(snum)

    # Open JLD2 file and unpack data corresonding to scenario snum
    jldfile = jldopen("13_pipelines_150_scenarios.jld2","r")
    junction_data = jldfile["junction_data$snum"]
    pipeline_data = jldfile["pipeline_data"]
    compressor_data = jldfile["compressor_data"]
    close(jldfile)

    # Setup topology dictionaries
    junction_map_in = Dict()    # Pipeline into each junction
    junction_map_out = Dict()   # Pipelines out of each junction
    pipe_map = Dict()           # Junction connected to each pipeline
    compressor_map = Dict()
    jmap = Dict()

    #create junction models.
    junctions = []
    for (i,j_data) in junction_data
        # Create a junction OptiGraph from j_data
        jmodel = create_junction_model(j_data,nt)
        # Create empty list with key corresponding to the OptiGraph
        # These will be filled with either a compressor or pipeline later
        junction_map_in[jmodel] = []
        junction_map_out[jmodel] = []
        # Add the junction OptiGraph to jmap and junctions
        jmap[i] = jmodel
        push!(junctions,jmodel)
    end

    pipelines = []
    compressors = []
    for (i,pdata) in pipeline_data
        #create pipeline Optigraph
        pipe_model = create_pipeline_model(pdata,nt,nx)
        # Create map of which junction is connected to which pipeline
        j_from = jmap[pdata[:from_node]]
        j_to = jmap[pdata[:to_node]]
        pipe_map[pipe_model] = [j_from,j_to]
        push!(junction_map_out[j_from], pipe_model)
        push!(junction_map_in[j_to], pipe_model)
        push!(pipelines,pipe_model)
    end

    for (i,cdata) in compressor_data
        #create compressor OptiGraph
        comp_model = create_compressor_model(cdata,nt)
        # Create map of which junctions are connected to which compressors
        j_from = jmap[cdata[:from_node]]
        j_to = jmap[cdata[:to_node]]
        compressor_map[comp_model] = [j_from,j_to]
        push!(junction_map_out[j_from],comp_model)
        push!(junction_map_in[j_to],comp_model)
        push!(compressors,comp_model)
    end

    # Create a general OptiGraph
    # This OptiGraph will be made of subgraphs of all the junctions, pipelines, and compressors
    gas_network = OptiGraph()

    # Add pipeline subgraphs to gas_network OptiGraph
    for pipe in pipelines
        add_subgraph!(gas_network,pipe)
    end

    # Add compressors to gas_network OptiGraph
    for compressor in compressors
        add_subgraph!(gas_network,compressor)
    end

    # Add junctions to gas_network OptiGraph
    for junction in junctions
        add_subgraph!(gas_network,junction)
    end

    # Create linking constraints between the compressors and the junctions on either side of them
    for compressor in compressors
        junction_from,junction_to = compressor_map[compressor]
        nodes = compressor[:time_nodes]
        @linkconstraint(gas_network,[t = 1:nt],nodes[t][:psuction] == junction_from[:time_nodes][t][:pressure])
        @linkconstraint(gas_network,[t = 1:nt],nodes[t][:pdischarge]  == junction_to[:time_nodes][t][:pressure])
    end

    # Create linking constraints and definitions for connections between the pipelines and junctions on either side of them 
    for junction in junctions
        devices_in = junction_map_in[junction]
        devices_out = junction_map_out[junction]

        if length(devices_in) > 0
            flow_in = [sum(devices_in[i][:fout][t] for i = 1:length(devices_in)) for t = 1:nt]
        else
            flow_in = zeros(nt)
        end
        if length(devices_out) > 0
            flow_out = [sum(devices_out[i][:fin][t] for i = 1:length(devices_out)) for t = 1:nt]
        else
            flow_out = zeros(nt)
        end
        total_supplied = [junction[:time_nodes][t][:total_supplied] for t = 1:nt]
        total_delivered = [junction[:time_nodes][t][:total_delivered] for t = 1:nt]

        @linkconstraint(gas_network,[t = 1:nt], flow_in[t] - flow_out[t] + total_supplied[t] - total_delivered[t] == 0)
    end

    # Fix demands within the scenario
    for (i,j_data) in junction_data
        jmodel = jmap[i]
        dvalues = j_data[:demand_values]
        n_demands = length(dvalues)
        nodes = jmodel[:time_nodes]
        for (t,node) in enumerate(nodes)
            @constraint(node,[d = 1:n_demands],node[:fdemand][d] == dvalues[d][t])
        end
    end

    # Return the gas network OptiGraph 
    return gas_network

end