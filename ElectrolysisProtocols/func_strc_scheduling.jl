#--------------------------------------------
# Define a function that schedules a structured operation of multiple devices over a single day. 
#--------------------------------------------

# Packages
using JuMP
using Gurobi
using DataFrames
using XLSX
using CSV

function strc_scheduling(date, N_stack, C_stack, protocols, minmax = false, stack_plot = true)
    """
        date: a specific operating date (1,2,3,4,5,6 or 7: July 1-7 2024)
        N_stack: number of electrolyzer stacks
        C_stack: maximum capacity of each electrolyzer in MW, identically applied to all stacks
        protocols: "v1" or "v2" or "v3" (details of each set of protocols are in data.jl)
        minmax: whether to ignore protocols or not, true or false
        stack_plot: whether to plot power consumption profile of individual stacks, true or false
    """
    #--------------------------------------------
    # Data and parameters
    #--------------------------------------------    
    global date
    # Load market data and protocols
    include("data.jl")

    if protocols == "v1"
        prot = prot_L_8 
    elseif protocols == "v2"
        prot = prot_diff_L
    elseif protocols == "v3"
        prot = prot_L_4
    end

    if length(prot) != N_stack
        error("Reassign protocols to $(N_stack) electrolyzers")
    end

 #   max_load_prot = [maximum(prot[k]) for k in K] # maximum load of protocol in %
    K = 1:N_stack 
    L_prot = [length(prot[k]) for k in K] # slot legnth of protocol 
    S_prot = [length(T) ÷ L_prot[k] for k in K] # number of slots over the entire time horizon

    # New time scale 
    N = Vector{Any}(undef, N_stack) 
    S = Vector{Any}(undef, N_stack) 
    lambda_M_new = Vector{Any}(undef, N_stack) 

    for k in K
        N[k] = 1:L_prot[k] # a set of inner-slot time steps
        S[k] = 1:S_prot[k] # a set of number of time slots
        lambda_M_new = reshape(lambda_M, L_prot[k],S_prot[k]) # reshape lambda_M to match the new time scale
    end

    # Electrolyzer parameters
    C_E = N_stack * C_stack # total system capacity in MW

    P_max = Vector{Any}(undef, N_stack)
    P_min = Vector{Any}(undef, N_stack)
    lambda_start = Vector{Any}(undef, N_stack)

    for k in K
#        P_max[k] = C_stack * max_load_prot[k] # maximum power load in MW
        P_max[k] = C_stack * max_load # maximum power load in MW
        P_min[k] = C_stack * min_load # minimum power load in MW
        lambda_start[k] = C_stack * start_cost # cold start-up cost in $
    end

    # Hydrogen storage
    C_S_out = C_E * eta_full_load # storage output flow capacity in kg/h
    C_S = C_E * eta_full_load * 24 * 1.0  # storage capacity in kg
    
    # Hydrogen demand
    D_daily = C_E * eta_full_load * 24 * 0.75 # daily hydrogen demand in kg

    #--------------------------------------------
    # Optimization
    #--------------------------------------------     

    model = Model(Gurobi.Optimizer)

    # Variables
    @variables model begin
        e[k in K, S[k], N[k]] >= 0     # power consumption by each device in MW
        h[k in K, S[k], N[k]] >= 0     # hydrogen production by each device in kg/h
        u_on[k in K, S[k]], Bin        # indicator for the usage of each device K at each time slot in S
        u_off[k in K, S[k]], Bin       # indicator for the non-usage of each device K at each time slot in S
        u_start[k in K, S[k]], Bin     # indicator for the start-up of each device K at each time slot in S
        u_min[k in K, S[k]], Bin       # indicator for the min usage of each device K at each time slot in S
        u_max[k in K, S[k]], Bin       # indicator for the max usage of each device K at each time slot in S
        u_prot[k in K, S[k]], Bin      # indicator for the protocol usage of each device K at each time slot in S

        e_Tdim[K, T] >= 0              # power consumption by each device defined in T scale in MW
        h_Tdim[K, T] >= 0              # hydrogen production from each device defined in T scale in kg/h

        e_tot[T] >= 0                  # power consumption by all devices in total in MW
        h_tot[T] >= 0                  # hydrogen production by all devices in total in kg/h
        h_d[T] >= 0      # hydrogen production delivered directly to demand in kg/h
        d[T] >= 0        # hydrogen sold in kg/h
        s_in[T] >= 0     # hydrogen injected into storage in kg/h
        s_out[T] >= 0    # hydrogen released from storage in kg/h
        soc[T] >= 0      # state of charge of storage in kg/h
    end

    # Objective: Minimize electricity cost
    @objective(model, Min,
        sum(e_tot[t] * lambda_M[t] for t=T) + 
        sum(sum(u_start[k,σ] * lambda_start[k] for σ=S[k]) for k=K)
    ) 
    
    # Constraints

    # Electrolyzers
    @constraint(model, [k=K, σ=S[k]], 
        u_on[k,σ] + u_off[k,σ] == 1) # exclusive states
    @constraint(model, [k=K, σ=S[k]],
        u_start[k,σ] >= (σ > 1 ? u_on[k,σ] - u_on[k,σ-1] : 0)) # cold start-up
    if minmax == false
        @constraint(model, [k=K, σ=S[k]], 
            u_max[k,σ] + u_min[k,σ] + u_prot[k,σ] == u_on[k,σ]) # exclusive operating modes with protocol
        @constraint(model, [k=K, σ=S[k], τ=N[k]], 
            e[k,σ,τ] == u_max[k,σ] * P_max[k] + u_min[k,σ] * P_min[k] + u_prot[k,σ] * prot[k][τ]) # power consumption in MW
    elseif minmax == true
        @constraint(model, [k=K, σ=S[k]], 
            u_max[k,σ] + u_min[k,σ] == u_on[k,σ]) # exclusive operating modes without protocol
        @constraint(model, [k=K, σ=S[k], τ=N[k]], 
            e[k,σ,τ] == u_max[k,σ] * P_max[k] + u_min[k,σ] * P_min[k]) # power consumption in MW
    end
    @constraint(model, [k=K, σ=S[k], τ=N[k]], 
        h[k,σ,τ] == A * e[k,σ,τ] + B * u_on[k,σ]) # hydrogen production in kg/h

    @constraint(model, [k=K, σ=S[k], τ=N[k]],
        e_Tdim[k, length(N[k]) * (σ - 1) + τ] == e[k,σ,τ]) # time scale mapping 
    @constraint(model, [k=K, σ=S[k], τ=N[k]],
        h_Tdim[k, length(N[k]) * (σ - 1) + τ] == h[k,σ,τ]) # time scale mapping
    
    @constraint(model, [t=T],
        e_tot[t] == sum(e_Tdim[k,t] for k=K)) # total power consumption in MW
    @constraint(model, [t=T],
        h_tot[t] == sum(h_Tdim[k,t] for k=K)) # total hydrogen production in kg/h

    # Hydrogen storage and demand
    @constraint(model, [t=T],
        h_tot[t] == h_d[t] + s_in[t]) # hydrogen balance in kg/h
    @constraint(model, [t=T],
        d[t] == h_d[t] + s_out[t]) # hydrogen balance in kg/h
    @constraint(model, [t=T],
        soc[t] == (t > 1 ? soc[t-1] : 0) + s_in[t] * time_res_rt - s_out[t] * time_res_rt) # storage balance in kg
    @constraint(model, [t=T],
        soc[t] <= C_S) # hydrogen storage fill in kg 
    @constraint(model, [t=T],
        s_out[t] <= C_S_out) # maximum storage output flow in kg/h
    @constraint(model,
        sum(d[t] * time_res_rt for t=T) >= D_daily) # daily hydrogen demand in kg

    # Solve
    start = time()
    set_silent(model)
    optimize!(model)

    #--------------------------------------------
    # Results
    #--------------------------------------------
    # Report results
    println("-------------------------------------");
    if termination_status(model) == MOI.OPTIMAL
        println(objective_value(model))
        println("\n\nRESULTS:")
        println("Objective cost         = $(round(objective_value(model), digits=digs)) \$")
        println("Minimum H2 demand      = $(round(D_daily, digits=digs)) kg")
        println("Amount of H2 demand    = $(round(sum(value.(d[t]) * time_res_rt for t=T), digits=digs)) kg")
        println("Max SOC                = $(maximum(value.(soc)))")
        println("Total computing time   = $(time() - start) sec")
        println("Nodes Explored         = $(node_count(model))")
        else
        println("No solution")
    end
    println("\n--------------------------------------");

    # Profile
    p_prof = [[value.(e_Tdim)[k,t] for t in T] for k in K] # power consumption by electrolyzer
    total_p_prof = sum(p_prof[k] for k in K) # total power consumption
    s_prof = [value.(soc)[t] for t in T] # hydrogen storage

    # Optimal cost
    optimal_cost = round(objective_value(model), digits=digs) # optimal cost in $

   # Plot 1. RTM and total power consumption profile
    plot(T * time_res_rt, lambda_M/10,
        label = "RTM [ X \$10/MWh]",
        ylabel = "RTM price [ X \$10/MWh]",
        line = 4,
        color = :slategray4,
        yaxis = :right
        )
    plot!(T * time_res_rt, total_p_prof,
        label = "Power",
        ylabel = "Power consumption [MW]",
        line = 4,
        color = :blue1
        )
    display(plot!(
        xlabel = "Time [h]",
        framestyle = :box,
        size = (600,300)
        ))

    # Plot 2. RTM and hydrogen storage profile
    plot(T * time_res_rt, lambda_M/100,
        label = "RTM [ X \$100/MWh]",
        ylabel = "RTM price [ X \$100/MWh]",
        line = 4,
        color = :slategray4,
        yaxis = :right
        )
    plot!(T * time_res_rt, s_prof / C_S,
        label = "H2",
        ylabel = "Hydrogen SOC [kg/kg]",
        line = 4,
        color = :green
        )
    display(plot!(
        xlabel = "Time [h]",
        framestyle = :box,
        size = (600,300)
        ))

    # Plot 3. Power consumption profiles of individual stacks
    if stack_plot == true
        for k in K
            plot(T * time_res_rt, p_prof[k,:],
            label = "Power",
            line = 4,
            color = :black
            )
        display(plot!(
            xlabel = "Time [h]",
            framestyle = :box,
            size = (800, 100),
            legend = :false
            ))
        end
    end

    return T, lambda_M, total_p_prof, s_prof, p_prof, optimal_cost
end
