#--------------------------------------------
# Define a function that schedules a structured operation of multiple devices over a single day. 
#--------------------------------------------

# Packages
using JuMP
using Gurobi
using DataFrames
using XLSX
using CSV

function flex_scheduling(date, N_stack, C_stack, R, baseline = false, stack_plot = true)
    """
        date: a specific operating date (1,2,3,4,5,6 or 7: July 1-7 2024)
        N_stack: number of electrolyzer stacks
        C_stack: maximum capacity of each electrolyzer in MW, identically applied to all stacks
        R: ramping rate between 0 to 1, identically applied to all stacks
        baseline: whether to use a baseline operation, true or false   
        stack_plot: whether to plot power consumption profile of individual stack, true or false
    """
    #--------------------------------------------
    # Data and parameters
    #--------------------------------------------    
    global date
    # Load market data and electrolyzer parameters
    include("data.jl")

    # Additional electrolyzer parameters
    K = 1:N_stack
    C_E = N_stack * C_stack # total systm capacity in MW

    P_max = Vector{Any}(undef, N_stack)
    P_min = Vector{Any}(undef, N_stack)
    lambda_start = Vector{Any}(undef, N_stack)

    for k in K
        P_max[k] = C_stack * max_load # maximum power load in MW
        P_min[k] = C_stack * min_load # minimum power load in MW
        lambda_start[k] = C_stack * start_cost # cold start-up cost in $
    end

    # Hydrogen storage
    C_S_out = C_E * eta_full_load # storage output flow capacity in kg/h
    C_S = C_E * eta_full_load * 24 * 1.0 # storage capacity in kg
    
    # Hydrogen demand
    D_daily = C_E * eta_full_load * 24 * 0.75 # daily hydrogen demand in kg

    #--------------------------------------------
    # Optimization
    #--------------------------------------------   

    model = Model(Gurobi.Optimizer)

    # Variables
    @variables model begin
        e[K, T] >= 0        # power consumption by each device in MW
        h[K, T] >= 0        # hydrogen production by each device in kg/h
        z_on[K,T], Bin      # electrolyzer state ON
        z_off[K,T], Bin     # electrolyzer state OFF
        z_start[K,T], Bin   # electrolyzer state START-UP
        
        e_tot[T] >= 0       # power consumption by all devices in total in MW
        h_tot[T] >= 0       # hydrogen production from all deivces in total in kg/h
        h_d[T] >= 0         # hydrogen production delivered directly to demand in kg/h
        d[T] >= 0           # hydrogen sold in kg/h
        s_in[T] >= 0        # hydrogen injected into storage in kg/h
        s_out[T] >= 0       # hydrogen released from storage in kg/h
        soc[T] >= 0         # state of charge of storage in kg/h
    end

    # Objective: Minimize electricity cost
    @objective(model, Min, 
        sum(e_tot[t] * lambda_M[t,1] for t=T) + 
        sum(z_start[k,t] * lambda_start[k] for t=T for k=K)
    )
    
    # Constraints

    # Electrolyzers
    @constraint(model, [k=K, t=T],
        1 == z_on[k,t] + z_off[k,t]) # exclusive states   
    @constraint(model, [k=K, t=T],
        z_start[k,t] >= (t > 1 ? z_on[k,t] - z_on[k,t-1] : 0)) # cold start-up
    if baseline == false   
        @constraint(model, [k=K, t=T],
            e[k,t] <= P_max[k] * z_on[k,t]) # maximum capacity in MW
        @constraint(model, [k=K, t=T],
            e[k,t] >= P_min[k] * z_on[k,t]) # minimum capacity in MW
        @constraint(model, [k=K, t=T[2:end]],
            e[k,t] - e[k,t-1] <= R * P_max[k]) # ramping 
        @constraint(model, [k=K, t=T[2:end]],
            e[k,t] - e[k,t-1] >= -R * P_max[k])  # ramping
    elseif baseline == true
        @constraint(model, [k=K, t=T],
            z_on[k,t] == 1) # baseline: always on
        @constraint(model, [k=K, t=T],
            e[k,t] == P_max[k]) # baseline: full load
    end
    @constraint(model, [k=K, t=T],
        h[k,t] == (A * e[k,t] + B * z_on[k,t])) # hydrogen production in kg/h

    
    @constraint(model, [t=T],
        e_tot[t] == sum(e[k,t] for k=K)) # total power consumption in MW
    @constraint(model, [t=T],
        h_tot[t] == sum(h[k,t] for k=K)) # total hydrogen production in kg/h

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
    p_prof = [[value.(e)[k,t] for t in T] for k in K] # power consumption by electrolyzer
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