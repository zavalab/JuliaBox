# Data Consolidation for the Different Formulations and Methods and compute the SSEs

# Load in all of the required packages
using FileIO, JLD, PyPlot, PyCall

# Load in all of the relevant data files
# Experimental Data
data_synthetic = load("./data/synthetic_data.jld2")
M = data_synthetic["Data"][:Monospecies]    # Monospecies experiment data
P = data_synthetic["Data"][:Pairwise]       # Pairwise experiment data     

# Model Results
# Extract the data for the discrete formulation
modeldata = load("./data/discrete_results.jld2")
model_output_d = modeldata["Data"]
# Extract the data for the continuous Formulations
modeldata = load("./data/continuous_results.jld2")
model_output_C = modeldata["Data"]

# Experimental Results
# Define system parameters          
n_m = length(M)        # Number of mono-species experiments and species
n_p = length(P)        # Number of pairwise experiments
n_p_s = length(P[1])   # Number of experiments within each pairwise experiment (There are more than one due to dilutions) 

# Include the concentration function
function get_concentration(d::Dict, t)
    ts = collect(keys(d))
    closest_index = findmin(abs.(ts.-t))[2]
    nearest_t = ts[closest_index]
    return d[nearest_t]
end

# Restructure the experimental data for plotting
exp_data = Dict()
# Monospecies experiments
for i = 1:n_m
    time = M[i][:supports]
    x_data = zeros(size(time))
    for j = 1:length(time)
        x_data[j] = M[i][:data][time[j]]
    end
    dict1 = Dict(:time => time, :data => x_data)
    exp_data[M[i][:species], M[i][:species]] = dict1
end
# Pairwise experiments
for i = 1:n_p
    for k = P[i][1][:species]
        t_1 = P[i][1][:supports]
        t_2 = t_1[end] .+ P[i][2][:supports]
        t_3 = t_2[end] .+ P[i][3][:supports]
        time = [t_1; t_2; t_3]
        x_1 = ones(size(t_1)); x_2 = ones(size(t_2)); x_3 = ones(size(t_3))
        for j = 1:length(t_1)
            x_1[j] = get_concentration(P[i][1][:data][k], t_1[j])
        end
        for j = 1:length(t_2)
            x_2[j] = get_concentration(P[i][2][:data][k], t_2[j] .- t_1[end])
        end
        for j = 1:length(t_3)
            x_3[j] = get_concentration(P[i][3][:data][k], t_3[j] .- t_2[end])
        end
        x_data = [x_1; x_2; x_3]
        dict1 = Dict(:time => time, :data => x_data)
        if k == P[i][1][:species][1]
            n = P[i][1][:species][2]
        else
            n = P[i][1][:species][1]
        end
        exp_data[k, n] = dict1
    end
end

# Restructure the real data
exp_real_data = Dict()
M = data_synthetic["Data"][:Monospecies_Real]    # Monospecies experiment data
P = data_synthetic["Data"][:Pairwise_Real]       # Pairwise experiment data     

# Monospecies experiments
for i = 1:n_m
    time = M[i][:supports]
    x_data = zeros(size(time))
    for j = 1:length(time)
        x_data[j] = M[i][:data][time[j]]
    end
    dict1 = Dict(:time => time, :data => x_data)
    exp_real_data[M[i][:species], M[i][:species]] = dict1
end
# Pairwise experiments
for i = 1:n_p
    for k = P[i][1][:species]
        t_1 = P[i][1][:supports]
        t_2 = t_1[end] .+ P[i][2][:supports]
        t_3 = t_2[end] .+ P[i][3][:supports]
        time = [t_1; t_2; t_3]
        x_1 = ones(size(t_1)); x_2 = ones(size(t_2)); x_3 = ones(size(t_3))
        for j = 1:length(t_1)
            x_1[j] = get_concentration(P[i][1][:data][k], t_1[j])
        end
        for j = 1:length(t_2)
            x_2[j] = get_concentration(P[i][2][:data][k], t_2[j] .- t_1[end])
        end
        for j = 1:length(t_3)
            x_3[j] = get_concentration(P[i][3][:data][k], t_3[j] .- t_2[end])
        end
        x_data = [x_1; x_2; x_3]
        dict1 = Dict(:time => time, :data => x_data)
        if k == P[i][1][:species][1]
            n = P[i][1][:species][2]
        else
            n = P[i][1][:species][1]
        end
        exp_real_data[k, n] = dict1
    end
end

# Store all of the data in one dictionary
alldata = Dict(:Discrete => model_output_d, :Continuous => model_output_C, 
                :Experimental => exp_data, :Experimental_Real => exp_real_data)

# Define the actual parameters
μ_actual = [0.245 0.246 0.584 0.237 0.478 0.457 0.598 0.402 0.219 0.502 0.232 0.156]
α_actual = [-0.912 0.453 0 0 0 0.137 0 0.692 0.961 0 0 1.343;                                      
            -0.306 -0.829 0 -0.56 0 -0.657 0 -1.098 0 -0.241 0.042 0;
            -0.228 -0.261 -0.88 -0.324 -0.632 -.584 -0.754 -0.124 0.231 -0.151 -0.176 -0.061;
            -0.529 -0.671 0 -0.622 0 0 0 -1.077 -0.402 -0.771 -0.433 0;
            -0.215 -0.278 -0.921 -0.265 -0.734 -0.556 -0.819 -0.114 -0.099 -0.465 -0.203 -0.020;
            -0.129 -0.168 -0.55 -0.202 -0.515 -0.66 -0.755 -0.049 0.759 0.039 -0.026 -0.030;
            -0.272 -0.274 -0.816 -0.303 -0.624 -0.642 -0.907 -0.169 -0.069 -0.632 -0.2 0;
            0.176 -0.448 3.378 -0.904 1.757 1.303 2.272 -2.442 -0.768 0.024 0.176 -0.139;
            -0.231 -1.121 -0.782 -0.405 -0.209 -0.638 -0.703 -0.153 -1.038 -0.507 0 -0.168;
            -0.348 0.312 0.066 0.265 -0.507 -0.047 -0.086 0 0.448 -1.454 -0.151 1.083;
            -0.9 0 0 -0.977 0 -0.106 0 -.408 1.015 -2.157 -1.254 0;
            -0.546 0 0 -0.817 0 0 -0.738 0 0 -0.439 0 -1.271]

# Extract the parameters for the different parameter estimation formulations
data = alldata["Data"]
SSE_α = Dict()
SSE_μ = Dict()
for i in keys(data)
    if i == :Discrete
        SSE_α[0] = sum((data[i][:Parameters][:α][n, m] - α_actual[n, m]) ^ 2 for n in 1:12 for m in 1:12)            # Index 0 references the discrete formulation
        SSE_μ[0] = sum((data[i][:Parameters][:μ][n] - μ_actual[n]) ^ 2 for n in 1:12)
    elseif i == :Continuous
        for j = 2:2:6
            SSE_α[j] = sum((data[i][j][:Parameters][:α][n, m] - α_actual[n, m]) ^ 2 for n in 1:12 for m in 1:12)     # Index 0 references the discrete formulation
            SSE_μ[j] = sum((data[i][j][:Parameters][:μ][n] - μ_actual[n]) ^ 2 for n in 1:12)
        end
    end
end

# Print the results
for i = 0:2:6
    print("SSE_μ_$i:")
    display(SSE_μ[i])
end
for i = 0:2:6
    print("SSE_α_$i:")
    display(SSE_α[i])
end
