# Parameter Estimation for a Dynamic, Biological System Using a Sigmoidal Functional Form for Each Experiment

# Load in all required packages
using FileIO, JLD2, PyPlot, LsqFit

# Load in the data file
datasets = load("./data/synthetic_data.jld2")
M = datasets["Data"][:Monospecies]     # Monospecies data
P = datasets["Data"][:Pairwise]        # Pairwise data

# Define the system parameters
n_m = 12                # Number of monospecies experiments
n_p = 66                # Number of pairwise experiments
n_p_s = 3               # Number of sub-experiments
p_n = 6                 # Number of parameters for each experiment

# Create a storage array for the monospecies parameters 
param_store_m = Dict()

# Perform parameter estimations for each pairwise experiment
for i = 1:n_m
    # Extract the data for the corresponding time points
    t_data = M[i][:supports]
    x_1_data = zeros(size(t_data))
    for k = 1:length(x_1_data)
        x_1_data[k] = M[i][:data][t_data[k]]
    end

    # Define the empirical model form, initial guess, and fit the model
    m(t, p) = p[1] ./(p[2] .+ p[3] .* exp.(p[4] .* (t .- p[5])))
    if i == 1
        p_0 = [1.0; 1.0; 1.0; 1.0; 1.0]
    else 
        p_0 = param_store_m[i-1]
    end 
    fit = curve_fit(m, t_data, x_1_data, p_0)
    
    # Extract and store the parameter values
    param = fit.param
    param_store_m[i] = param
    display(param)
    # Plot the results to assess the fit
    pygui(true)
    subplot(4, 3, i) 
    plot(t_data, x_1_data, ".")
    plot(t_data[1]:1:t_data[end], m(t_data[1]:1:t_data[end], param))
end

# Create a storage array for the pairwise parameters
param_store_p = Dict()

# Perform parameter estimations for each pairwise experiment
for i = 1:n_p
    for j = 1:n_p_s
        # Extract the data for the corresponding time point and species 
        t_data = P[i][j][:supports]
        x_1_data = zeros(size(P[i][j][:supports]))
        x_2_data = zeros(size(P[i][j][:supports]))
        for k = 1:length(P[i][j][:supports])
            x_1_data[k] = P[i][j][:data][P[i][j][:species][1]][P[i][j][:supports][k]]
            x_2_data[k] = P[i][j][:data][P[i][j][:species][2]][P[i][j][:supports][k]]
        end

        # Define the empirical model form, initial guess, and fit the model
        m(t, p) = p[1] ./(p[2] .+ p[3] .* exp.(p[4] .* (t .- p[5])))
        if i == 1 && j ==1
            p_0_1 = [1.0; 1.0; 1.0; 1.0; 1.0]
            p_0_2 = [1.0; 1.0; 1.0; 1.0; 1.0]
        elseif i == 28 || i == 29
            p_0_1 = [1.0; 1.0; 1.0; 1.0; 1.0]
            p_0_2 = [1.0; 1.0; 1.0; 1.0; 1.0]
        else 
            p_0_1 = param_store_p[1, 1, P[1][j][:species][1]]
            p_0_2 = param_store_p[1, 1, P[1][j][:species][2]]
        end

        fit1 = curve_fit(m, t_data, x_1_data, p_0_1)
        fit2 = curve_fit(m, t_data, x_2_data, p_0_2)

        # Extract and store the parameter values
        param1 = fit1.param
        param2 = fit2.param
        param_store_p[i, j, P[i][j][:species][1]] = param1
        param_store_p[i, j, P[i][j][:species][2]] = param2

        # Plot the results to assess the fit
        # pygui(true)
        # n = (i - 1) * n_p_s  + j
        # subplot(18, 11, n) 
        # plot(t_data, x_1_data, ".")
        # plot(t_data[1]:1:t_data[end], m(t_data[1]:1:t_data[end], param1))
        # plot(t_data, x_2_data, ".")
        # plot(t_data[1]:1:t_data[end], m(t_data[1]:1:t_data[end], param2))
    end
end

# Store all of the results
param_store = Dict(:Monospecies => param_store_m, :Pairwise => param_store_p)

# Save the parameter values to a data file 
# save("./data/sigmoidal_params.jld2", "Data", param_store)