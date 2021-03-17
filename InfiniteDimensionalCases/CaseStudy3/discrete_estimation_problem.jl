# Parameter Estimation for a Dynamic, Biological System Using InfiniteOpt.jl with Discrete Synthetic Data

# Load in all required packages
using FileIO, JLD, JLD2, JuMP, InfiniteOpt, Ipopt

# Load in the data files
data_synthetic = load("./data/synthetic_data.jld2")
M = data_synthetic["Data"][:Monospecies]    # Monospecies experiment data
P = data_synthetic["Data"][:Pairwise]       # Pairwise experiment data

# Include the parameter estimation function
include("estimation_model.jl")

# Parameter estimation arguments
args = Dict(:Mono_Data => M,                                    # Assign monospecies data
            :Pair_Data => P,                                    # Assign pairwise data
            :Formulation => :Discrete,                          # Formulation (:Discrete or :Continuous)
            :Derivative_Method => OrthogonalCollocation(2))

# Solve the model
results = Param_Est_Function(args)

# Save the results of the discrete formulation
save("discrete_results.jld2", "Data", results)