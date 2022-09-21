# Import packages required for running model
using Plasmo
using MadNLP
using MadNLPGraph
using MadNLPHSL
using JLD2
using DelimitedFiles


#######################################################################################

#Model Data
horizon = 24*3600        # The time horizon in seconds
nt= 24                   # Number of time points
nx = 10                  # Number of space points per pipeline
dt = horizon / (nt - 1)  # Time delta
NS = 10                  # Number of scenarios
S = 1:NS

# Load in the function to build the pipelines
include("modelfunctions_stoch.jl")
# Load in the data
include("load_data_stoch.jl")

#######################################################################################

# Make overall optigraph
gas_network_stoch = OptiGraph()

# Add a master node to contain first stage variables
@optinode(gas_network_stoch,master)

# Set the first stage variable to be compressor power
# Variable is indexed by compressor number time points
@variable(master,power[1:11,1:nt])

# Create array for subgraphs; each subgraph will contain a scenario
gas_networks = Array{Any,1}(undef,NS)


for l in 1:NS
    # Generate a subgraph containing one scenario
    gas_networks[l] = make_stoch_subgraph(l)
    # Add that scenario subgraph to the overall OptiGraph
    add_subgraph!(gas_network_stoch,gas_networks[l])
end

# Create a set of all subgraphs within the overall OptiGraph
stoch_subs = getsubgraphs(gas_network_stoch)

# Create array for the set of subgraphs within each scenario
scenario_subs       = Array{Any,1}(undef,NS)

for j in S
    # Get the set of all subgraphs within each scenario
    scenario_subs[j] = getsubgraphs(stoch_subs[j])
end

for i in 1:11
    for j in 1:24
        for k in S
            # Link the first stage variables
            # Indices are k = scenario number, i = compressor number, and j = time points
            # Within scenario_subs, there are 13 subgraphs for pipelines, and then the compressor subgraphs
            # so the compressor subgraph call below is offset by 13
            @linkconstraint(gas_network_stoch, scenario_subs[k][i+13][:time_nodes][j][:power] == master[:power][i,j])
        end
    end
end

# Aggregate the model into a single OptiGraph with no subgraphs
# This aggregates every subgraph into a node, so gas_node below has NS nodes + a master node
gas_node,dict = aggregate(gas_network_stoch,0)

gas_node_agg, dict_agg = aggregate(gas_network_stoch)


#######################################################################################

# Solve the problem with MadNLP.optimize!
# This code uses Ma57 as the linear solver, but this solver requires additional setup
# If you just want to run this example problem without Ma57, you can use linear_solver=MadNLPUmfpack instead in line 82 below
MadNLP.optimize!(gas_node; linear_solver=MadNLPMa57)
#MadNLP.optimize!(gas_node; linear_solver=MadNLPSchur, schur_custom_partition=true, schur_subproblem_solver=MadNLPMa57)

#######################################################################################

# Save the optimized values of gas_node to a CSV

nodes      = all_nodes(gas_node)
all_values = Array{Any,2}(undef, (11376,11))

# Set the first column to be the variable name
all_values[:,1] = all_variables(nodes[1])

# Set columns 2 - NS + 1 to be the value of the variable in each scenario
for i in 1:NS
    vars = all_variables(nodes[i])
    for j in 1:11376
        all_values[j, i + 1] = value(vars[j])
    end
end

master_node_vars   = all_variables(nodes[NS + 1])
master_node_values = value.(master_node_vars)

master_node = Array{Any,2}(undef, (264,2))
master_node[:,1] = master_node_vars
master_node[:,2] = master_node_values

DelimitedFiles.writedlm("10_scenario_results.csv", all_values, ',')
DelimitedFiles.writedlm("10_scenario_master_node.csv", master_node, ',')

#######################################################################################
# Below is optional code for plotting a single scenario of the gas network. This was used to create the visualizations of the graph network
# The locations of the nodes and edges from the generated graph was used to produce clearer visualizaitons

using PlasmoPlots

# PlasmoPlots.jl is available at: https://github.com/jalving/PlasmoPlots.jl.git


network = make_stoch_subgraph(1)

plt1 = layout_plot(network, node_labels=false,  markersize=5,linewidth=.1,linealpha=2.0,
subgraph_colors=true, layout_options=Dict(:tol => .05, :C=>100, :K =>.005,  :iterations=> 400 ),
plt_options=Dict(:size=>(500,500),:legend=>false, :framestyle=>:box, :grid => false, :axis => nothing));

display(plt1)
