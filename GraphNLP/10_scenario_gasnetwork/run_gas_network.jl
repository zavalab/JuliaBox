# Import packages required for running model
using Plasmo
using MadNLP
using MadNLPGraph
using MadNLPHSL
using JLD2

#######################################################################################

#Model Data
horizon = 24*3600        # The time horizon in seconds
nt= 24                   # Number of time points
nx = 10                  # Number of space points per pipeline
dt = horizon / (nt - 1)  # Time delta
NS = 5                   # Number of scenarios
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
#MadNLP.optimize!(gas_node; linear_solver=MadNLPMa57)
MadNLP.optimize!(gas_node; linear_solver=MadNLPMa57)
#MadNLP.optimize!(gas_node; linear_solver=MadNLPSchur, schur_custom_partition=true, schur_subproblem_solver=MadNLPMa57)