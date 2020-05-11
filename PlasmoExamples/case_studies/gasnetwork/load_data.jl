using JLD2

jldfile = jldopen("13_pipelines.jld2","r")
junction_data = jldfile["junction_data"]
pipeline_data = jldfile["pipeline_data"]
compressor_data = jldfile["compressor_data"]
close(jldfile)

#setup topology dictionaries
junction_map_in = Dict()    #pipeline into each junction
junction_map_out = Dict()   #pipelines out of each junction
pipe_map = Dict()           #junction connected to each pipeline
compressor_map = Dict()
jmap = Dict()

#create junction models.
junctions = []
for (i,j_data) in junction_data
    jmodel = create_junction_model(j_data,nt)
    junction_map_in[jmodel] = []
    junction_map_out[jmodel] = []
    jmap[i] = jmodel
    push!(junctions,jmodel)
end

pipelines = []
compressors = []
for (i,pdata) in pipeline_data
    #create pipe models
    pipe_model = create_pipeline_model(pdata,nt,nx)
    j_from = jmap[pdata[:from_node]]
    j_to = jmap[pdata[:to_node]]
    pipe_map[pipe_model] = [j_from,j_to]
    push!(junction_map_out[j_from], pipe_model)
    push!(junction_map_in[j_to], pipe_model)
    push!(pipelines,pipe_model)
end

for (i,cdata) in compressor_data
    #create compressor models
    comp_model = create_compressor_model(cdata,nt)
    j_from = jmap[cdata[:from_node]]
    j_to = jmap[cdata[:to_node]]
    compressor_map[comp_model] = [j_from,j_to]
    push!(junction_map_out[j_from],comp_model)
    push!(junction_map_in[j_to],comp_model)
    push!(compressors,comp_model)
end
