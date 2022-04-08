using DelimitedFiles

# Load in the variables from each scenario graph (11376 x 10 array)
all_values    = DelimitedFiles.readdlm("10_scenario_results.csv", ',')[:,2:11]
# Load in the variables from the master node (length 264 array)
master_values = DelimitedFiles.readdlm("10_scenario_master_node.csv", ',')[:,3]

# Note that the order of the compressor and pipelines is not sequential because the OptiGraph was built from dictionary data
# The lists below are the order that the compressors and pipelines show up in the solution
compressor_order = [5,7,8,1,4,6,11,2,10,9,3]
pipe_order       = [5,7,12,8,1,4,6,13,11,2,10,9,3]

# px_vals is a 4-D array; The index order is scenario, pipeline number, time point, and spacial discretization point
# This array contains the pressure at each discretized point of the pipelines; this pressure is used in the linepack calculation
# Note that the order of pipelines is given by "pipe_order" above. So px_vals[:,1,:,:] is the 5th pipeline since pipe_order[1] = 5
px_vals = Array{Any,4}(undef, (10,13,24,10))
for i in 1:10
    for j in 1:13
        for l in 1:10
            for k in 1:24
                # Indexing is challenging because of the number of variables in this problem. 
                # See the CSV for examples of how the variables are sorted
                index = (j-1)*240*3 + 1 + (k-1) * 3 + (l-1)*24*3
                px_vals[i,j,k,l] = all_values[index,i]
            end
        end
    end
end

# Define the linepack, which is the average of the first nine discretized point
# This only uses first nine points to be consistent with the linepack constraint of the model
# Note that the linepack is an average pressure (in bars)
# Linepack can be converted to mass by multiplying by the volume of the pipeline (D = .92 m, L is pipeline dependent; see JLD2 file)
# and dividing by the speed of sound squared (see make_JLD2_stoch_150.jl for this value, which is given as "nu2" in units of m^2/s^2)
linepack = Array{Any,3}(undef, (10,13,24))

for i in 1:13
    for j in 1:10
        for k in 1:24
            linepack[j,i,k] = sum(px_vals[j,i,k,l] for l in 1:9)/9
        end
    end
end

# Define an array for the compressor power; rows represent the time from 1-24 and columns represent the compressor
# Note that the compressor order is given by "compressor_order" above
# Units of compressor power are kW
comp_power = Array{Any,2}(undef, (24,11))

for j in 1:24
    for i in 1:11
        comp_power[j,i] = master_values[(j-1)*11 + i]
    end
end


# Define arrays for the amount of gas delivered and the demand values
# Rows represent the scenario number and columns represent the time point
# Units are SCMx10^4/hr for each time period; this can be converted to mass by using 
# the density value given in "make_JLD2_stoch_150.jl"
deliver_values = Array{Any,2}(undef, (10,24))
demand_values  = Array{Any,2}(undef, (10,24))

for j in 1:24
    for i in 1:10
        index = 11114 + (j-1)*4
        deliver_values[i,j] = all_values[index, i]
        demand_values[i,j]  = all_values[(index+1),i]
    end
end
