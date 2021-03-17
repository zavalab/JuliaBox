using JuMP, InfiniteOpt, Ipopt, JLD

function Param_Est_Function(args::Dict)
    # Load in the data files that are independent of the formulation
    data_sets = load("./data/source_data.jld")
    A = data_sets["Data"][:A] # Interaction parameter indices, bounds, and start values
    U = data_sets["Data"][:U] # Growth parameter, indices, bounds, and start values

    # Assign the monospecies and pairwise data to their respective variables
    M = args[:Mono_Data] # Mono-species experiment data
    P = args[:Pair_Data] # Pairwise experiment data

    # Set the number of supports
    if args[:Formulation] == :Continuous
        Param_Mono = args[:Empirical_Parameters][:Monospecies]
        Param_Pair = args[:Empirical_Parameters][:Pairwise]
    end
    
    # Define system parameters          
    n_m = length(M) # Number of mono-species experiments and species
    n_p = length(P) # Number of pairwise experiments
    n_p_s = length(P[1]) # Number of experiments within each pairwise experiment (There are more than one due to dilutions) 
    number_supports = args[:Number_Supports] # Number of supports

    # Define the infinite model
    model = InfiniteModel(Ipopt.Optimizer)
    set_optimizer_attribute(model, "max_iter", 5000)

    # Define the parameters to be estimated (finite optimization variables)
    # Growth rate parameters U and interaction parameters A
    @hold_variable(model, U[:μ_lb][i] ≤ μ[i ∈ 1:n_m] ≤ U[:μ_ub][i], 
                   start = U[:μ_start][i])
    @hold_variable(model, A[:α_lb][i, j] ≤ α[i ∈ 1:n_m, j ∈ 1:n_m] ≤ A[:α_ub][i, j], 
                   start = A[:α_start][i, j])

    # Define the infinite time parameters for each experiment
    if args[:Formulation] == :Continuous
        # Mono-species experiments (only define one infinite parameter since each experiment has identical supports)
        @infinite_parameter(model, t_m ∈ [M[1][:t_0], M[1][:t_f]], num_supports = number_supports, 
                            derivative_method = args[:Derivative_Method])
        # Pairwise experiments (define an infinite parameter for each indiviual experiment since they have different supports)
        @infinite_parameter(model, t_p[i ∈ 1:n_p, j ∈ 1:n_p_s] ∈ [P[i][j][:t_0], P[i][j][:t_f]], num_supports = number_supports, 
                            independent = true, derivative_method = args[:Derivative_Method])
    elseif args[:Formulation] == :Discrete
        # Mono-species experiments (only define one infinite parameter since each experiment has identical supports)
        @infinite_parameter(model, t_m ∈ [M[1][:t_0], M[1][:t_f]], supports = M[1][:supports], 
                            derivative_method = args[:Derivative_Method])
        # Pairwise experiments (define an infinite parameter for each indiviual experiment since they have different supports)
        @infinite_parameter(model, t_p[i ∈ 1:n_p, j ∈ 1:n_p_s] ∈ [P[i][j][:t_0], P[i][j][:t_f]], supports = P[i][j][:supports], 
                            independent = true, derivative_method = args[:Derivative_Method])
    end

    # Define the infinite variables 
    # Mono-species concentrations x_m[experiment, species](t)
    @infinite_variable(model, 0 ≤ x_m[i ∈ 1:n_m, j ∈ M[i][:species]](t_m), 
                       start = M[i][:init_cond][1])
    # Mono-species dummy variables α[species] * x_m[experiment, species](t)
    @infinite_variable(model, α_x_m[i=1:n_m, j=M[i][:species]](t_m), start = 0)
    # Pairwise concentrations x_p[experiment, sub-experiment, species](t)
    @infinite_variable(model, 0 ≤ x_p[i ∈ 1:n_p, j ∈ 1:n_p_s, k ∈ P[i][j][:species]](t_p[i, j]), 
                       start = P[i][j][:init_cond][k])
    # Pairwise dummy variables α[species1, species2] * x_p[experiment, sub-experiment, species1, species2](t)
    @infinite_variable(model, α_x_p[i ∈ 1:n_p, j ∈ 1:n_p_s, k ∈ P[i][j][:species], n=P[i][j][:species]](t_p[i, j]), 
                       start = 0)

    # Define constraints for each experiment
    # Mono-species experiments
    @constraint(model, [i ∈ 1:n_m, j ∈ M[i][:species]], α_x_m[i, j] == x_m[i, j] * α[j, j])
    @constraint(model, [i ∈ 1:n_m, j ∈ M[i][:species]], 
                ∂(x_m[i, j], t_m) == (μ[j] + α_x_m[i, j]) * x_m[i, j])
    # Pairwise experiments
    @constraint(model, [i ∈ 1:n_p, j ∈ 1:n_p_s, k ∈ P[i][j][:species], n ∈ P[i][j][:species]], 
                α_x_p[i, j, k, n] == α[n, k] * x_p[i, j, k])
    @constraint(model, [i ∈ 1:n_p, j ∈ 1:n_p_s, k ∈ P[i][j][:species][1], n ∈ P[i][j][:species][2]], 
                ∂(x_p[i, j, k], t_p[i, j]) == (μ[k] + α_x_p[i, j, k, k] + α_x_p[i, j, n, k]) * x_p[i, j, k])  
    @constraint(model, [i ∈ 1:n_p, j ∈ 1:n_p_s, k ∈ P[i][j][:species][2], n ∈ P[i][j][:species][1]], 
                ∂(x_p[i, j, k], t_p[i, j]) == (μ[k] + α_x_p[i, j, k, k] + α_x_p[i, j, n, k]) * x_p[i, j, k])

    # Define parameter functions for the experimental data, which is dependent on the formulation type
    if args[:Formulation] == :Continuous
        # Mono-species experiments
        x_m_d = Dict()
        for i=1:n_m
            p = Param_Mono[i]
            function get_concentration(t)
                return p[1] ./(p[2] .+ p[3] .* exp.(p[4] .* (t .- p[5])))
            end
            x_m_d[i, M[i][:species]] = parameter_function(get_concentration, t_m)
        end
        # Pairwise experiments
        x_p_d = Dict()
        for i=1:n_p
            for j=1:n_p_s
                for k = P[i][j][:species]
                    p = Param_Pair[i, j, k]
                    function get_concentration(t)
                        return p[1] ./(p[2] .+ p[3] .* exp.(p[4] .* (t .- p[5])))
                    end
                    x_p_d[i, j, k] = parameter_function(get_concentration, t_p[i, j])
                end
            end
        end
    else
        # Mono-species experiments
        x_m_d = Dict()
        for i=1:n_m
            d = M[i][:data]
            function get_concentration(t)
                ts = collect(keys(d))
                closest_index = findmin(abs.(ts.-t))[2]
                nearest_t = ts[closest_index]
                return d[nearest_t]
            end
            x_m_d[i, M[i][:species]] = parameter_function(get_concentration, t_m)
        end
        # Pairwise experiments
        x_p_d = Dict()
        for i=1:n_p
            for j=1:n_p_s
                for k = P[i][j][:species]
                    d = P[i][j][:data][k]
                    function get_concentration(t)
                        ts = collect(keys(d))
                        closest_index = findmin(abs.(ts.-t))[2]
                        nearest_t = ts[closest_index]
                        return d[nearest_t]
                    end
                    x_p_d[i, j, k] = parameter_function(get_concentration, t_p[i, j])
                end
            end
        end
    end

    # Calculate the sum of the squared residuals for the objective function
    # Monospecies experiments
    res_m = sum(support_sum((x_m_d[i, j] - x_m[i, j])^2, t_m) for i in 1:n_m for j in M[i][:species]) 
    # Pairwise experiments
    res_p = sum(support_sum((x_p_d[i, j, k] - x_p[i, j, k])^2 , t_p[i, j]) for i in 1:n_p for j in 1:n_p_s for k in P[i][j][:species][1:end]) 

    # Define the objective function
    @objective(model, Min, res_m + res_p)

    # Optimize the model
    optimize!(model)

    # Extract and save the results for later plotting
    model_output = Dict()
    # Monospecies experiments
    for i = 1:n_m
        time = supports(t_m)
        x_hat = JuMP.value.(x_m[i, M[i][:species]])
        dict1 = Dict(:time => time, :data => x_hat)
        model_output[M[i][:species], M[i][:species]] = dict1
    end
    # Pairwise experiments
    for i = 1:n_p
        for k = P[i][1][:species]
            t_1 = supports(t_p[i, 1])
            t_2 = t_1[end] .+ supports(t_p[i, 2])
            t_3 = t_2[end] .+ supports(t_p[i, 3])
            time = [t_1; t_2; t_3]
            x_hat = [JuMP.value.(x_p[i, 1, k]); JuMP.value.(x_p[i, 2, k]); JuMP.value.(x_p[i, 3, k])]
            dict1 = Dict(:time => time, :data => x_hat)
            if k == P[i][1][:species][1]
                j = P[i][1][:species][2]
            else
                j = P[i][1][:species][1]
            end
            model_output[k, j] = dict1
        end
    end
    # Extract the parameters and initial conditions for comparison between formulation type and derivative methods
    μ_result = JuMP.value.(μ)
    α_result = JuMP.value.(α)
    init_cond = Dict()
    for i = 1:n_p
        for j = 1:n_p_s
            for k = P[i][j][:species]
                init_cond[i, j, k] = JuMP.value.(x_p[i, j, k])[1]
            end
        end
    end
    param = Dict(:μ => μ_result, :α => α_result, :init_cond => init_cond)
    # Save the solution time
    solution_time = solve_time(model)
    # Save the model status
    status = raw_status(model)

    # Store all of the data
    results = Dict(:Model_Data => model_output, :Parameters => param, 
                   :Solve_Time => solution_time, :Status => status)

    # Return the results
    return results
end
