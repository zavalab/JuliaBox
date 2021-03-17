# Generating Synthetic Data for Parameter Estimation Using Simulations 

# Load in all of the relevant packages
using DifferentialEquations, DiffEqParamEstim, Optim, Distributions, FileIO, JLD, 
JLD2

# Literature Values for Growth and Interaction Parameters (gLV model) 
# Source: Venturelli OS, Carr AC, Fisher G, Hsu RH, Lau R, Bowen BP, et al. 
# "Deciphering microbial interactions in synthetic human gut microbiome communities. 
# Molecular Systems Biology". Appendix Figure S28
# Growth Rate
μ = [0.245 0.246 0.584 0.237 0.478 0.457 0.598 0.402 0.219 0.502 0.232 0.156]
# Interaction Parameters
α = [-0.912 0.453 0 0 0 0.137 0 0.692 0.961 0 0 1.343;                                      
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

# Load in the original data file from "Scalable Nonlinear Programming Framework 
# for Parameter Estimation in Dynamic Biological System Models"
# This is used to retrieve the experimental conditions
# We still need to generate synthetic data to obtain more data points
datasets = load("./data/source_data.jld")
M = datasets["Data"][:M]             # Monospecies experiment data
P = datasets["Data"][:P]             # Pairwise experiment data

# Choose a value for the standard deviation of the synthetic data
σ = 0.01

# Define system parameters
n_m = 12
n_p = 66
n_p_s = 3

# Create a storage dictionary for the monospecies synthetic data 
M_syn = Dict()
M_syn_real = Dict()

for i = 1:n_m
    # Extract the parameters for the experiment number
    μ_1 = μ[M[i][:species]]  
    α_11 = α[M[i][:species], M[i][:species]]

    # Create a function for the set of differential equations 
    function derivative_function(dx_dt, x, p, t)
        dx_dt[1] = (μ_1 .+ α_11 * x[1]) .* x[1]
    end

    # Define the time horizon, initial condition, and solve the problem
    tspan = (M[i][:t_0], M[i][:t_f])
    init_cond = [M[i][:init_cond]]
    prob = ODEProblem(derivative_function, init_cond, tspan)
    sol = DifferentialEquations.solve(prob, saveat = 1)

    # Extract the solution 
    x_1 = zeros(size(sol.t))
    x_1_real = zeros(size(sol.t))
    for i = 1:length(sol.t)
        x_1[i] = max(0, sol.u[i][1] + rand(Normal(0, σ)))
        x_1_real[i] = sol.u[i][1]
    end
    x_data = Dict(sol.t[j] => x_1[j] for j in 1:length(sol.t))
    x_data_real = Dict(sol.t[j] => x_1_real[j] for j in 1:length(sol.t))
    M_syn[i] = Dict(:species => M[i][:species], :init_cond => init_cond, 
                    :t_0 => M[i][:t_0], :t_f => M[i][:t_f] .+ 0.0001, 
                    :supports => sol.t, :data => x_data)
    M_syn_real[i] = Dict(:species => M[i][:species], :init_cond => init_cond, 
                         :t_0 => M[i][:t_0], :t_f => M[i][:t_f] .+ 0.0001, 
                         :supports => sol.t, :data => x_data_real)
end

# Create a storage dictionary for the pairwise synthetic data
P_syn = Dict()
P_syn_real = Dict()
P_syn_dil = Array{Dict{Any, Any}}(undef, 66, 3)
P_syn_dil_real = Array{Dict{Any, Any}}(undef, 66, 3)

# Generate synthetic data for each pairwise experiment
for i = 1:n_p
   for j = 1:n_p_s
        # Extract the parameters for the experiment number
        μ_1 = μ[P[i][j][:species][1]]  
        μ_2 = μ[P[i][j][:species][2]]
        α_11 = α[P[i][j][:species][1], P[i][j][:species][1]]
        α_21 = α[P[i][j][:species][2], P[i][j][:species][1]] 
        α_22 = α[P[i][j][:species][2], P[i][j][:species][2]] 
        α_12 = α[P[i][j][:species][1], P[i][j][:species][2]]

        # Create a function for the set of differential equations 
        function derivative_function(dx_dt, x, p, t)
            dx_dt[1] = (μ_1 .+ α_11 * x[1] .+ α_12 * x[2]) .* x[1]
            dx_dt[2] = (μ_2 .+ α_21 * x[1] .+ α_22 * x[2]) .* x[2]
        end

        # Define the time horizon, initial condition, and solve the problem
        tspan = (P[i][j][:t_0], P[i][j][:t_f])
        init_cond = [P[i][j][:init_cond][P[i][j][:species][1]], 
                     P[i][j][:init_cond][P[i][j][:species][2]]]
        prob = ODEProblem(derivative_function, init_cond, tspan)
        sol = DifferentialEquations.solve(prob, AutoTsit5(Rosenbrock23()), 
                                          saveat = 2, maxiters = 1e5)

        # Extract the solution 
        x_1 = zeros(size(sol.t))
        x_2 = zeros(size(sol.t))
        x_1_real = zeros(size(sol.t))
        x_2_real = zeros(size(sol.t))
        for i = 1:length(sol.t)
            x_1[i] = max(0, sol.u[i][1] + rand(Normal(0, σ)))
            x_2[i] = max(0, sol.u[i][2] + rand(Normal(0, σ)))
            x_1_real[i] = sol.u[i][1]
            x_2_real[i] = sol.u[i][2]
        end

        # Create data dictionaries for the data and actual solution
        x_dict = Dict(P[i][j][:species][1] => x_1, P[i][j][:species][2] => x_2)
        x_dict_real = Dict(P[i][j][:species][1] => x_1_real, P[i][j][:species][2] => x_2_real)
        init_dict = Dict(P[i][j][:species][k] => P[i][j][:init_cond][P[i][j][:species][k]] for k in 1:2)
        data_dict = Dict(P[i][j][:species][k] => Dict(sol.t[n] => x_dict[P[i][j][:species][k]][n] for n in 1:length(sol.t)) for k in 1:2)
        data_dict_real = Dict(P[i][j][:species][k] => Dict(sol.t[n] => x_dict_real[P[i][j][:species][k]][n] for n in 1:length(sol.t)) for k in 1:2)
        P_syn_dil[i, j] = Dict(:species => P[i][j][:species], :t_0 => P[i][j][:t_0], 
                               :t_f => P[i][j][:t_f] .+ 0.0001, :supports => sol.t, 
                               :init_cond => init_dict, :data => data_dict)
        P_syn_dil_real[i, j] = Dict(:species => P[i][j][:species], :t_0 => P[i][j][:t_0], 
                                    :t_f => P[i][j][:t_f] .+ 0.0001, :supports => sol.t, 
                                    :init_cond => init_dict, :data => data_dict_real)
    end
    # Store the experimental data and the actual solution
    P_syn[i] = P_syn_dil[i, :]
    P_syn_real[i] = P_syn_dil_real[i, :]
end

# Store all data for both experiment types
data = Dict(:Monospecies => M_syn, :Pairwise => P_syn, 
            :Monospecies_Real => M_syn_real, :Pairwise_Real => P_syn_real)

# Save the synthetic data
# save("./data/synthetic_data.jld2", "Data", data)
