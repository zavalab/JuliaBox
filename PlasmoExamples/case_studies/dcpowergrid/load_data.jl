using LightGraphs, Random, LinearAlgebra, Printf,DelimitedFiles, SparseArrays, PowerModels

function get_data(path)
    data = PowerModels.parse_file(path)
    A = calc_susceptance_matrix(data)

    args = Dict()
    args[:N] = size(A.matrix,1)   #number of buses
    args[:y] = A.matrix;
    args[:del] = spzeros(args[:N],args[:N])
    args[:g] = Graph(args[:N])    #the lightgraph
    for e in values(data["branch"])
        i = A.bus_to_idx[e["f_bus"]]
        j = A.bus_to_idx[e["t_bus"]]
        LightGraphs.add_edge!(args[:g],i,j)
        args[:del][i,j] = e["angmin"]
        args[:del][j,i] = e["angmin"]
    end

    #Voltage angle limits
    args[:val] =-ones(args[:N])*Inf
    args[:vau] = ones(args[:N])*Inf
    args[:ref] =  A.bus_to_idx[reference_bus(data)["index"]]   #reference buses
    args[:val][args[:ref]] = 0; args[:vau][args[:ref]] = 0

    gens = bus_gen_lookup(data["gen"],data["bus"])
    args[:c1]  = [[1e6,-1e6] for i=1:args[:N]]
    args[:c2]  = [[.0,.0] for i=1:args[:N]]
    args[:sl]= [[.0,-1e2] for i=1:args[:N]]
    args[:su]= [[1e2,.0] for i=1:args[:N]]
    args[:ng] = 2*ones(Int64,args[:N])


    for i=1:args[:N]
        bus = A.idx_to_bus[i]
        for j=1:length(gens[bus])
            if length(gens[bus][j]["cost"])!=0
                push!(args[:c1][i],gens[bus][j]["cost"][1])
                push!(args[:c2][i],gens[bus][j]["cost"][2])
                # push!(args[:c2][i],gens[bus][j]["cost"][2])
                push!(args[:sl][i],gens[bus][j]["pmin"])
                push!(args[:su][i],gens[bus][j]["pmax"])
                args[:ng][i] += 1
            end
        end
    end

    args[:Ng]=sum(args[:ng])

    args[:sd] = zeros(args[:N])
    for v in values(data["load"])
        args[:sd][A.bus_to_idx[v["load_bus"]]] = v["pd"]
    end

    return args
end

path = (@__DIR__)*"/pglib_opf_case9241_pegase.m"
args = get_data(path)
