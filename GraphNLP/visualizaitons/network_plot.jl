using Plots
using DelimitedFiles
using Colors

cols = Colors.distinguishable_colors(15)

# Note that this code uses positions that were found through PlasmoPlots.jl: https://github.com/jalving/PlasmoPlots.jl.git

### Code for graphing 1 scenario of the stochastic gas network ###
######################################################################################

# Read in node and edge locations
nodes = readdlm("Graph_Nodes_network.csv",',')
edges = readdlm("Graph_Edges_network.csv",',')

# Define which nodes are pipelines, juncitons, or compressors
node_positions = nodes[:,2:3]
pipeline_nodes   = node_positions[1:3120,:]
junction_nodes   = node_positions[3385:3984,:]
compressor_nodes = node_positions[3121:3384,:]
# Define marker size
m_size = 5

# Define plot
plt_1scenario = Plots.scatter(axis=nothing, size=(700,700),legend=(.1,.87),framestyle=:box,legendfontsize=8)

# Add to plot
Plots.scatter!(pipeline_nodes[:,1], pipeline_nodes[:,2],markersize=m_size,color=:deepskyblue,markerstrokewidth=.3,markerstrokecolor=:black,label="Pipeline Nodes") #or deepskyblue
Plots.scatter!(junction_nodes[:,1], junction_nodes[:,2],markersize=m_size,color=:orange,markerstrokewidth=.3,markerstrokecolor=:black,label="Junction Nodes")
Plots.scatter!(compressor_nodes[:,1], compressor_nodes[:,2],markersize=m_size,color=:firebrick,markerstrokewidth=.3,markerstrokecolor=:black, label="Compressor Nodes")

# Add edges (link constraints) to plot
for i in 1:length(edges[:,1])
    n_from = Int(edges[i,1])
    n_to   = Int(edges[i,2])
    if i == 1
        Plots.plot!(plt_1scenario, [node_positions[n_from,1], node_positions[n_to,1]], [node_positions[n_from,2], node_positions[n_to,2]]; linewidth=1, linecolor=:blue,linealpha=.3, label="Linking Constraint")
    else
        Plots.plot!(plt_1scenario, [node_positions[n_from,1], node_positions[n_to,1]], [node_positions[n_from,2], node_positions[n_to,2]];label=false, linewidth=1, linecolor=:blue,linealpha=.3 )
    end
end




### Code for graphing 4 scenarios of the stochastic gas network and the master node ###
#######################################################################################################


# Read in node and edge locations
nodes = readdlm("Graph_Nodes_network.csv",',')
edges = readdlm("Graph_Edges_network.csv",',')

# Define node positions for each scenario. The center of the plot is at (15,15), so we will offset each scenario in different directions
node_positions = nodes[:,2:3]
node_positions2 = nodes[:,2:3]
node_positions3 = nodes[:,2:3]
node_positions4 = nodes[:,2:3]

# Offset scenarios
node_positions2[:,2] = node_positions[:,2] .+ 30
node_positions3[:,1] = node_positions[:,1] .+ 30
node_positions4      = node_positions .+ 30

# Define marker size
m_size = 3

# Define pipeline, junction and compressor nodes for each scenario
pipeline_nodes   = node_positions[1:3120,:]
junction_nodes   = node_positions[3385:3984,:]
compressor_nodes = node_positions[3121:3384,:]

pipeline_nodes2   = node_positions[1:3120,:]
junction_nodes2   = node_positions[3385:3984,:]
compressor_nodes2 = node_positions[3121:3384,:]

pipeline_nodes2[:,2] = pipeline_nodes[:,2] .+ 30
junction_nodes2[:,2] = junction_nodes[:,2] .+ 30
compressor_nodes2[:,2] = compressor_nodes[:,2] .+ 30

pipeline_nodes3   = node_positions[1:3120,:]
junction_nodes3   = node_positions[3385:3984,:]
compressor_nodes3 = node_positions[3121:3384,:]

pipeline_nodes3[:,1] = pipeline_nodes[:,1] .+ 30
junction_nodes3[:,1] = junction_nodes[:,1] .+ 30
compressor_nodes3[:,1] = compressor_nodes[:,1] .+ 30

pipeline_nodes4   = pipeline_nodes .+ 30
junction_nodes4   = junction_nodes .+ 30
compressor_nodes4 = compressor_nodes .+ 30


# Read in edge data
edges2 = readdlm("Graph_Edges_network.csv",',')
edges3 = readdlm("Graph_Edges_network.csv",',')
edges4 = readdlm("Graph_Edges_network.csv",',')

# Offset edge data
edges2[:,2] = edges[:,2] .+ 30
edges3[:,1] = edges[:,1] .+ 30
edges4      = edges .+ 30


# Define plot
plt_4scenarios = Plots.scatter(axis=nothing, size=(700,700),legend=(.1,.85),framestyle=:box,legendfontsize=8)#legendfont=("Helvetica",8))

# Add Pipelines and junctions for each scenario
Plots.scatter!(pipeline_nodes[:,1], pipeline_nodes[:,2],markersize=m_size,color=:deepskyblue,label=:none,markerstrokewidth=.3,markerstrokecolor=:black,) #or deepskyblue
Plots.scatter!(junction_nodes[:,1], junction_nodes[:,2],markersize=m_size,color=:orange,label=:none,markerstrokewidth=.3,markerstrokecolor=:black)

Plots.scatter!(pipeline_nodes2[:,1], pipeline_nodes2[:,2], markersize=m_size,color=:deepskyblue,label=:none,markerstrokewidth=.3,markerstrokecolor=:black)
Plots.scatter!(junction_nodes2[:,1], junction_nodes2[:,2], markersize=m_size, color=:orange,label=:none,markerstrokewidth=.3,markerstrokecolor=:black)

Plots.scatter!(pipeline_nodes3[:,1], pipeline_nodes3[:,2], markersize=m_size,color=:deepskyblue,label=:none,markerstrokewidth=.3,markerstrokecolor=:black)
Plots.scatter!(junction_nodes3[:,1], junction_nodes3[:,2], markersize=m_size, color=:orange,label=:none,markerstrokewidth=.3,markerstrokecolor=:black)

Plots.scatter!(pipeline_nodes4[:,1], pipeline_nodes4[:,2], markersize=m_size,color=:deepskyblue,label=:none,markerstrokewidth=.3,markerstrokecolor=:black)
Plots.scatter!(junction_nodes4[:,1], junction_nodes4[:,2], markersize=m_size, color=:orange,label=:none,markerstrokewidth=.3,markerstrokecolor=:black)

# Add compressors for each scenario; this requires the compressors to be on top of the graph so they are more visible
Plots.scatter!(compressor_nodes[:,1], compressor_nodes[:,2],markersize=m_size,color=:firebrick,label=:none,markerstrokewidth=.3,markerstrokecolor=:black)
Plots.scatter!(compressor_nodes2[:,1], compressor_nodes2[:,2], markersize=m_size, color=:firebrick,label=:none,markerstrokewidth=.3,markerstrokecolor=:black)
Plots.scatter!(compressor_nodes3[:,1], compressor_nodes3[:,2], markersize=m_size, color=:firebrick,label=:none,markerstrokewidth=.3,markerstrokecolor=:black)
Plots.scatter!(compressor_nodes4[:,1], compressor_nodes4[:,2], markersize=m_size, color=:firebrick,label=:none,markerstrokewidth=.3,markerstrokecolor=:black)

# Add edges (linking constraints) to master node
for i in 1:length(compressor_nodes[:,1])
    Plots.plot!(plt_4scenarios, [15,compressor_nodes[i,1]], [15,compressor_nodes[i,2]]; linewidth=.1,linecolor=:blue,linealpha=.2,label=false)
    Plots.plot!(plt_4scenarios, [15, compressor_nodes2[i,1]], [15, compressor_nodes2[i,2]], linewidth=.1, linecolor=:blue, linealpha=.2, label=false)
    Plots.plot!(plt_4scenarios, [15, compressor_nodes3[i,1]], [15, compressor_nodes3[i,2]], linewidth=.1, linecolor=:blue, linealpha=.2, label=false)
    Plots.plot!(plt_4scenarios, [15, compressor_nodes4[i,1]], [15, compressor_nodes4[i,2]], linewidth=.1, linecolor=:blue, linealpha=.2, label=false)
end

# Add master node; there are also other nodes defined under the master node to make labeling of the legend easier
Plots.scatter!([15],[15], markersize=m_size, markercolor=:deepskyblue, markerstrokewidth=.3,markerstrokecolor=:black, label="Pipeline Nodes")
Plots.scatter!([15],[15], markersize=m_size, markercolor=:orange, markerstrokecolor=:black,markerstrokewidth=.3, label="Junction Nodes")
Plots.scatter!([15],[15], markersize=m_size, markercolor=:firebrick, markerstrokecolor=:black,markerstrokewidth=.3, label="Compressor Nodes")
Plots.scatter!([15],[15], markersize=m_size, markercolor=:grey, markerstrokecolor=:black,markerstrokewidth=.3, label="Master Node")
Plots.scatter!([15],[15], markersize=10, markercolor=:grey,label=:false)


# Add linking constraints
# This code takes some time, so an if statement is added to infrom the user how far the for loop has gone
for i in 1:length(edges[:,1])

    n_from = Int(edges[i,1])
    n_to   = Int(edges[i,2])

    if i == 1
        Plots.plot!(plt_4scenarios, [node_positions[n_from,1], node_positions[n_to,1]], [node_positions[n_from,2], node_positions[n_to,2]]; linewidth=.05, linecolor=:blue,linealpha=.15, label="Linking Constraint")
    else
        Plots.plot!(plt_4scenarios, [node_positions[n_from,1], node_positions[n_to,1]], [node_positions[n_from,2], node_positions[n_to,2]];label=false, linewidth=.05, linecolor=:blue,linealpha=.15 )
    end

    n_from2 = Int(edges[i,1])
    n_to2   = Int(edges[i,2])
    Plots.plot!(plt_4scenarios, [node_positions2[n_from2,1], node_positions2[n_to2,1]], [node_positions2[n_from2,2], node_positions2[n_to2,2]];label=false, linewidth=.05, linecolor=:blue,linealpha=.15 )
    
    n_from3 = Int(edges[i,1])
    n_to3   = Int(edges[i,2])
    Plots.plot!(plt_4scenarios, [node_positions3[n_from3,1], node_positions3[n_to3,1]], [node_positions3[n_from3,2], node_positions3[n_to3,2]];label=false, linewidth=.05, linecolor=:blue,linealpha=.15 )
    
    n_from4 = Int(edges[i,1])
    n_to4   = Int(edges[i,2])
    Plots.plot!(plt_4scenarios, [node_positions4[n_from4,1], node_positions4[n_to4,1]], [node_positions4[n_from4,2], node_positions4[n_to4,2]];label=false, linewidth=.05, linecolor=:blue,linealpha=.15 )


    if i%1000 == 0
        println("Done with iteration ", i, " out of ", length(edges[:,1]))
    end
end



### Code for graphing the stochastic PID problem in both time and scenario ###
######################################################################################################################


# Read in node and edge locations
nodes = readdlm("Graph_Nodes_PID.csv",',')
edges = readdlm("Graph_Edges_PID.csv",',')
node_positions = nodes[:,2:3]

# Choose colors
cols = Colors.distinguishable_colors(6)
if cols[1] == Colors.parse(Colorant,"black")
    cols[1] = Colors.parse(Colorant,"grey")
end
cols

# Define nodes corresponding to each scenario
scenario_1 = node_positions[2:101,:]
scenario_2 = node_positions[102:201,:]
scenario_3 = node_positions[202:301,:]
scenario_4 = node_positions[302:401,:]
scenario_5 = node_positions[402:501,:]

# Define nodes corresponding to each time
time_1 = vcat(node_positions[1:26,:], node_positions[102:126,:],node_positions[202:226,:], node_positions[302:326,:], node_positions[402:426,:])
time_2 = vcat(node_positions[27:51,:], node_positions[127:151,:],node_positions[227:251,:], node_positions[327:351,:], node_positions[427:451,:])
time_3 = vcat(node_positions[52:76,:], node_positions[152:176,:],node_positions[252:276,:], node_positions[352:376,:], node_positions[452:476,:])
time_4 = vcat(node_positions[77:101,:], node_positions[177:201,:], node_positions[277:301,:], node_positions[377:401,:], node_positions[477:501,:])

# Define plot for the scenario-partitioned graph
plt_PID_scenario = Plots.scatter(scenario_1[:,1], scenario_1[:,2],markercolor=cols[2],axis=nothing, legend=false, size=(500,500),framestyle=:box)
Plots.scatter!(scenario_2[:,1], scenario_2[:,2],markercolor=cols[3])
Plots.scatter!(scenario_3[:,1], scenario_3[:,2], markercolor=cols[4])
Plots.scatter!(scenario_4[:,1], scenario_4[:,2], markercolor=cols[5])
Plots.scatter!(scenario_5[:,1], scenario_5[:,2], markercolor=cols[6])
Plots.scatter!([node_positions[1,1]],[node_positions[1,2]], markercolor=cols[1])

# Plot edges (linking constraints)
for i in 1:length(edges[:,1])
    n_from = Int(edges[i,1])
    n_to   = Int(edges[i,2])
    Plots.plot!(plt_PID_scenario, [node_positions[n_from,1], node_positions[n_to,1]], [node_positions[n_from,2], node_positions[n_to,2]];label=false, linewidth=1, linecolor=:blue,linealpha=200 )
end

# Define plot for the time-partitioned graph
plt_PID_time = Plots.scatter(time_1[:,1], time_1[:,2],markercolor=cols[2],axis=nothing,legend=false, size=(500,500),framestyle=:box)
Plots.scatter!(time_2[:,1], time_2[:,2],markercolor=cols[3])
Plots.scatter!(time_3[:,1], time_3[:,2],markercolor=cols[4])
Plots.scatter!(time_4[:,1], time_4[:,2],markercolor=cols[5])

# Plot edges (linking constraints)
for i in 1:length(edges[:,1])
    n_from = Int(edges[i,1])
    n_to   = Int(edges[i,2])
    Plots.plot!(plt_PID_time, [node_positions[n_from,1], node_positions[n_to,1]], [node_positions[n_from,2], node_positions[n_to,2]];label=false, linewidth=1, linecolor=:blue,linealpha=200 )
end

