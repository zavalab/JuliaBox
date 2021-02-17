# Parameter Estimation for a Dynamic, Biological System Using the Continuous Formulation with InfiniteOpt.jl 

# Load in all required packages
using FileIO, JLD, JLD2, JuMP, InfiniteOpt, Ipopt

# Load in the data files
data_synthetic = load("./data/synthetic_data.jld2")
parameters = load("./data/sigmoidal_params.jld2")
M = data_synthetic["Data"][:Monospecies]    # Monospecies experiment data
P = data_synthetic["Data"][:Pairwise]       # Pairwise experiment data
empirical_parameters = parameters["Data"]   # Parameters for the empirical function

# Include the parameter estimation function
include("estimation_model.jl")

# Create a dictionary for storing the results
results = Dict()

# Solve the model with orthogonal collocation using different number of nodes
for i = 2:2:6
    # Parameter estimation arguments
    args = Dict(:Mono_Data => M,                                    # Assign monospecies data
                :Pair_Data => P,                                    # Assign pairwise data
                :Formulation => :Continuous,                        # Formulation (:Discrete or :Continuous)
                :Number_Supports => 15,                             # Number of supports (only important for continuous formulation)
                :Empirical_Parameters => empirical_parameters,      # Empirical parameters (only important for continuous formulation)
                :Derivative_Method => OrthogonalCollocation(i))
    results[i] = Param_Est_Function(args)
end

# Save the results of the continuous formulation
save("continuous_results.jld2", "Data", results)