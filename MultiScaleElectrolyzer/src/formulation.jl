using JuMP, Gurobi, DataFrames, XLSX, Plots
using CSV
include("electrolyzer_struct.jl")

"""
Function to run the model as a 3-state problem (on/off/sb)
"""
function run_3st_opt(θ::Electrolyzer, spp_array; gap=0.01)
    model = Model(Gurobi.Optimizer)
    ϕ = θ.ϕ # nominal capacity
    α_max = θ.α_max
    α_init = α_max
    α_min = θ.α_min
    δ_on = θ.δ_on
    δ_start = θ.δ_start
    ρ = θ.i
    B = 9.66
    ρ_sb = θ.ρ_sb
    ΔT = 1
    n_years = Int.(ceil(length(spp_array)/8760))
    years = Int.(1:(n_years))

    spp_avg_arr = spp_array
    T_span = collect(1:length(spp_avg_arr))

    # this will contain only the points 8760, 2*8760, 3*8760 ...
    T_star = [t+1 for (i, t) in enumerate(T_span) if i % (8760) == 0] # 720 --> for once a month in DAM, 2880 --> once a month in 15-RTM
    
    # this is T_span \ {T_star U t1}
    T_subset = setdiff(T_span, T_star)
    T_subset = setdiff(T_subset, [1])

    # Create a list where each element is the range of t for that year
    year_indices = [((y-1)*8760 + 1):(y*8760) for y in years]
    year_indices[end] = ((n_years-1)*8760+1):(length(spp_avg_arr)) # in case spp_arr is not exactly 22*8760

    # new ------------
    @variables model begin
        e_tot[T_span] >= 0 # electrolyzer consumption 
        h[T_span] >= 0 # hydrogen production
        z_on[T_span], Bin # on electrolyzer
        z_sb[T_span], Bin # on electrolyzer
        z_off[T_span], Bin # on electrolyzer 
        z_start[T_span], Bin # relaxed start binary 
        z_replace[years], Bin
        A[T_span] >= 0
        w[T_span] >= 0 # bilinearization of on binary
        v[years] >= 0 # bilinearization of replacement binary
    end


    # Device model constraints
    @constraint(model, [t in T_span], e_tot[t] == ϕ * z_on[t] * ΔT + ϕ * z_sb[t] * ΔT * ρ_sb) 

    # Efficiency constraints
    @constraint(model, A[1] == α_init - δ_on * z_on[1] * ΔT - δ_start * z_start[1]) # IC
    @constraint(model, [t in T_subset], A[t] == A[t-1] - δ_on * z_on[t] * ΔT - δ_start * z_start[t])
    @constraint(model, [(y,t) in zip(years, T_star)], A[t] == A[t-1] - v[y] - δ_on * z_on[t] * ΔT - δ_start * z_start[t] + α_max * z_replace[y]) # replacements


    # Bilinear term linearization (v_t = A_{t-1} * z_replace_{t} ) # if replace = 1 --> A[t-1] cancels out, otherwise =0
    @constraint(model, [y in years], v[y] <= α_max * z_replace[y])
    @constraint(model, [y in years], v[y] >= α_min * z_replace[y])
    @constraint(model, [(y,t) in zip(years, T_star)], v[y] <= A[t-1] - α_min * (1 - z_replace[y]))
    @constraint(model, [(y,t) in zip(years, T_star)], v[y] >= A[t-1] - α_max * (1 - z_replace[y]))

    @constraint(model, [t in T_span], h[t] == ϕ * w[t] * ΔT + B * z_on[t] * ΔT)

    # Bilinear term linearization (ω_t = A_t * z_t)
    @constraint(model, [t in T_span], w[t] <= α_max * z_on[t])
    @constraint(model, [t in T_span], w[t] >= α_min * z_on[t])
    @constraint(model, [t in T_span], w[t] >= A[t] - α_max * (1 - z_on[t]))
    @constraint(model, [t in T_span], w[t] <= A[t] - α_min * (1 - z_on[t]))

    # 3 state logic
    @constraint(model, [t in T_span], z_on[t] + z_off[t] + z_sb[t] == 1)

    # meet a daily demand quota
    n_hours_per_day = 24
    total_days = Int.(floor(length(T_span)/n_hours_per_day))
    D_H = 750 #sigma
    daily_indices = [
        T_span[(d-1)*n_hours_per_day + 1 : d*n_hours_per_day]
        for d in 1:total_days
    ]

    # 3. Build the constraint by iterating over your new groups
    @constraint(model, daily_quota[i in daily_indices],
        sum(h[t] for t in i) >= D_H
    )
    # # Start up binary logic
    @constraint(model, [t=2:length(T_span)], z_start[t] >= z_on[t] - z_on[t-1] - z_sb[t-1] ) # binary constraints
    @constraint(model, [t=2:length(T_span)], z_start[t] <= z_on[t])
    @constraint(model, [t=2:length(T_span)], z_start[t] <= 1 - z_on[t-1] - z_sb[t-1])

    # Economic-specific parameters
    λ_H = θ.λ_H # hydrogen selling price [$/kg]
    λ_CAPEX_yearly = θ.λ_CAPEX_yearly
    λ_CAPEX = θ.λ_CAPEX_Plant
    λ_OPEX = θ.λ_OPEX
    λ_CAPEX_replace = θ.λ_CAPEX_Stack

    @expression(model, cap_exp, λ_CAPEX * ϕ) # [=] $/MW * MW
    @expression(model, op_exp[y in years], λ_OPEX)
    @expression(model, replace_exp[y in years], ϕ * λ_CAPEX_replace * z_replace[y])
    # @expression(model, h_profit, sum(h[t] * λ_H for t in T_span)) # [=] kg H2 * $/kg H2
    @expression(model, h_revenue[y in year_indices], sum(h[t] * λ_H for t in y))
    @expression(model, e_exp[y in year_indices], sum(e_tot[t] * spp_avg_arr[t] for t in y)) # [=] MWhr * $/MWhr
    # @expression(model, start_costs, sum(z_start[t] * start_cost for t in T_span))
    # @expression(model, profit, h_profit - elec_exp - cap_exp - op_exp - replace_exp) # [=] $
    
    # TODO: make sets for the indices of t belonging to each year
    @expression(model, NPV[(y, y_k) in zip(years, year_indices)], (1/(1+ρ)^(y-1))* (-replace_exp[y] - e_exp[y_k] - op_exp[y] + h_revenue[y_k]))
    # Maximize profit
    @objective(model, Min, -sum(NPV[y] for y in zip(years, year_indices)))
    @expression(model, value, sum(NPV[y] for y in zip(years, year_indices)) - cap_exp)
    set_optimizer_attribute(model, "MIPGap", gap)
    # set_optimizer_attribute(model, "IterationLimit", 500000)
    # set_optimizer_attribute(model, "TimeLimit", 500)

    @expression(model, factor[y in years], (1/(1+ρ)^(y-1)))
    @expression(model, LCOH_num, λ_CAPEX + sum(factor[y]*(e_exp[y_k] + op_exp[y] + replace_exp[y]) for (y, y_k) in zip(years, year_indices)))
    @expression(model, LCOH_den, sum(factor[y]*sum(h[t] for t in y_k) for (y, y_k) in zip(years, year_indices)))
    @expression(model, LCOH, LCOH_num/LCOH_den)


    optimize!(model)

    return model
end

function run_3st_opt_constant(θ::Electrolyzer, spp_array; gap=0.01)
    model = Model(Gurobi.Optimizer)
    ϕ = θ.ϕ # nominal capacity
    α_max = θ.α_max
    α_init = α_max
    α_min = θ.α_min
    δ_on = θ.δ_on
    δ_start = θ.δ_start
    ρ = θ.i
    B = 9.66
    ρ_sb = θ.ρ_sb
    ΔT = 1
    n_years = Int.(ceil(length(spp_array)/8760))
    years = Int.(1:(n_years))

    spp_avg_arr = spp_array
    T_span = collect(1:length(spp_avg_arr))

    # this will contain only the points 8760, 2*8760, 3*8760 ...
    T_star = [t+1 for (i, t) in enumerate(T_span) if i % (8760) == 0] # 720 --> for once a month in DAM, 2880 --> once a month in 15-RTM
    
    # this is T_span \ {T_star U t1}
    T_subset = setdiff(T_span, T_star)
    T_subset = setdiff(T_subset, [1])

    # Create a list where each element is the range of t for that year
    year_indices = [((y-1)*8760 + 1):(y*8760) for y in years]
    year_indices[end] = ((n_years-1)*8760+1):(length(spp_avg_arr)) # in case spp_arr is not exactly 22*8760

    # new ------------
    @variables model begin
        e_tot[T_span] >= 0 # electrolyzer consumption 
        h[T_span] >= 0 # hydrogen production
        z_on[T_span], Bin # on electrolyzer
        z_sb[T_span], Bin # on electrolyzer
        z_off[T_span], Bin # on electrolyzer 
        z_start[T_span], Bin # relaxed start binary 
        z_replace[years], Bin
        A[T_span] >= 0
        w[T_span] >= 0 # bilinearization of on binary
        v[years] >= 0 # bilinearization of replacement binary
    end


    # Device model constraints
    @constraint(model, [t in T_span], e_tot[t] == ϕ * z_on[t] * ΔT + ϕ * z_sb[t] * ΔT * ρ_sb) 

    # Efficiency constraints
    @constraint(model, A[1] == α_init - δ_on * z_on[1] * ΔT - δ_start * z_start[1]) # IC
    @constraint(model, [t in T_subset], A[t] == A[t-1] - δ_on * z_on[t] * ΔT - δ_start * z_start[t])
    @constraint(model, [(y,t) in zip(years, T_star)], A[t] == A[t-1] - v[y] - δ_on * z_on[t] * ΔT - δ_start * z_start[t] + α_max * z_replace[y]) # replacements


    # Bilinear term linearization (v_t = A_{t-1} * z_replace_{t} ) # if replace = 1 --> A[t-1] cancels out, otherwise =0
    @constraint(model, [y in years], v[y] <= α_max * z_replace[y])
    @constraint(model, [y in years], v[y] >= α_min * z_replace[y])
    @constraint(model, [(y,t) in zip(years, T_star)], v[y] <= A[t-1] - α_min * (1 - z_replace[y]))
    @constraint(model, [(y,t) in zip(years, T_star)], v[y] >= A[t-1] - α_max * (1 - z_replace[y]))

    @constraint(model, [t in T_span], h[t] == ϕ * w[t] * ΔT + B * z_on[t] * ΔT)

    # Bilinear term linearization (ω_t = A_t * z_t)
    @constraint(model, [t in T_span], w[t] <= α_max * z_on[t])
    @constraint(model, [t in T_span], w[t] >= α_min * z_on[t])
    @constraint(model, [t in T_span], w[t] >= A[t] - α_max * (1 - z_on[t]))
    @constraint(model, [t in T_span], w[t] <= A[t] - α_min * (1 - z_on[t]))

    # 3 state logic
    @constraint(model, [t in T_span], z_on[t] + z_off[t] + z_sb[t] == 1)
    @constraint(model, [t in T_span], z_on[t] == 1) # ALWAYS ON!!!!

    # meet a daily demand quota
    n_hours_per_day = 24
    total_days = Int.(floor(length(T_span)/n_hours_per_day))
    D_H = 750 #sigma
    daily_indices = [
        T_span[(d-1)*n_hours_per_day + 1 : d*n_hours_per_day]
        for d in 1:total_days
    ]

    # 3. Build the constraint by iterating over your new groups
    @constraint(model, daily_quota[i in daily_indices],
        sum(h[t] for t in i) >= D_H
    )
    # # Start up binary logic
    @constraint(model, [t=2:length(T_span)], z_start[t] >= z_on[t] - z_on[t-1] - z_sb[t-1] ) # binary constraints
    @constraint(model, [t=2:length(T_span)], z_start[t] <= z_on[t])
    @constraint(model, [t=2:length(T_span)], z_start[t] <= 1 - z_on[t-1] - z_sb[t-1])

    # Economic-specific parameters
    λ_H = θ.λ_H # hydrogen selling price [$/kg]
    λ_CAPEX_yearly = θ.λ_CAPEX_yearly
    λ_CAPEX = θ.λ_CAPEX_Plant
    λ_OPEX = θ.λ_OPEX
    λ_CAPEX_replace = θ.λ_CAPEX_Stack

    @expression(model, cap_exp, λ_CAPEX * ϕ) # [=] $/MW * MW
    @expression(model, op_exp[y in years], λ_OPEX)
    @expression(model, replace_exp[y in years], ϕ * λ_CAPEX_replace * z_replace[y])
    # @expression(model, h_profit, sum(h[t] * λ_H for t in T_span)) # [=] kg H2 * $/kg H2
    @expression(model, h_revenue[y in year_indices], sum(h[t] * λ_H for t in y))
    @expression(model, e_exp[y in year_indices], sum(e_tot[t] * spp_avg_arr[t] for t in y)) # [=] MWhr * $/MWhr
    # @expression(model, start_costs, sum(z_start[t] * start_cost for t in T_span))
    # @expression(model, profit, h_profit - elec_exp - cap_exp - op_exp - replace_exp) # [=] $
    
    # TODO: make sets for the indices of t belonging to each year
    @expression(model, NPV[(y, y_k) in zip(years, year_indices)], (1/(1+ρ)^(y-1))* (-replace_exp[y] - e_exp[y_k] - op_exp[y] + h_revenue[y_k]))
    
    @expression(model, factor[y in years], (1/(1+ρ)^(y-1)))
    @expression(model, LCOH_num, λ_CAPEX + sum(factor[y]*(e_exp[y_k] + op_exp[y] + replace_exp[y]) for (y, y_k) in zip(years, year_indices)))
    @expression(model, LCOH_den, sum(factor[y]*sum(h[t] for t in y_k) for (y, y_k) in zip(years, year_indices)))
    @expression(model, LCOH, LCOH_num/LCOH_den)

    # Maximize profit
    @objective(model, Min, -sum(NPV[y] for y in zip(years, year_indices)))
    @expression(model, value, sum(NPV[y] for y in zip(years, year_indices)) - cap_exp)
    set_optimizer_attribute(model, "MIPGap", gap)
    # set_optimizer_attribute(model, "IterationLimit", 500000)
    # set_optimizer_attribute(model, "TimeLimit", 500)
    optimize!(model)

    return model
end



"""
Function to run the model as a 2-state problem (on/off)
"""
function run_2st_opt(θ::Electrolyzer, spp_array; gap=0.01)
    model = Model(Gurobi.Optimizer)
    ϕ = θ.ϕ # nominal capacity
    α_max = θ.α_max
    α_init = α_max
    α_min = θ.α_min
    δ_on = θ.δ_on
    δ_start = θ.δ_start
    ρ = θ.i
    B = 9.66
    ρ_sb = θ.ρ_sb
    ΔT = 1
    n_years = Int.(ceil(length(spp_array)/8760))
    years = Int.(1:(n_years))

    spp_avg_arr = spp_array
    T_span = collect(1:length(spp_avg_arr))

    # this will contain only the points 8760, 2*8760, 3*8760 ...
    T_star = [t+1 for (i, t) in enumerate(T_span) if i % (8760) == 0] # 720 --> for once a month in DAM, 2880 --> once a month in 15-RTM
    
    # this is T_span \ {T_star U t1}
    T_subset = setdiff(T_span, T_star)
    T_subset = setdiff(T_subset, [1])

    # Create a list where each element is the range of t for that year
    year_indices = [((y-1)*8760 + 1):(y*8760) for y in years]
    year_indices[end] = ((n_years-1)*8760+1):(length(spp_avg_arr)) # in case spp_arr is not exactly 22*8760

    # new ------------
    @variables model begin
        e_tot[T_span] >= 0 # electrolyzer consumption 
        h[T_span] >= 0 # hydrogen production
        z_on[T_span], Bin # on electrolyzer
        z_off[T_span], Bin # off electrolyzer 
        z_start[T_span], Bin # start binary 
        z_replace[years], Bin
        A[T_span] >= 0
        w[T_span] >= 0 # bilinearization of on binary
        v[years] >= 0 # bilinearization of replacement binary
    end


    # Device model constraints
    @constraint(model, [t in T_span], e_tot[t] == ϕ * z_on[t] * ΔT) 

    # Efficiency constraints
    @constraint(model, A[1] == α_init - δ_on * z_on[1] * ΔT - δ_start * z_start[1]) # IC
    @constraint(model, [t in T_subset], A[t] == A[t-1] - δ_on * z_on[t] * ΔT - δ_start * z_start[t])
    @constraint(model, [(y,t) in zip(years, T_star)], A[t] == A[t-1] - v[y] - δ_on * z_on[t] * ΔT - δ_start * z_start[t] + α_max * z_replace[y]) # replacements


    # Bilinear term linearization (v_t = A_{t-1} * z_replace_{t} ) # if replace = 1 --> A[t-1] cancels out, otherwise =0
    @constraint(model, [y in years], v[y] <= α_max * z_replace[y])
    @constraint(model, [y in years], v[y] >= α_min * z_replace[y])
    @constraint(model, [(y,t) in zip(years, T_star)], v[y] <= A[t-1] - α_min * (1 - z_replace[y]))
    @constraint(model, [(y,t) in zip(years, T_star)], v[y] >= A[t-1] - α_max * (1 - z_replace[y]))

    @constraint(model, [t in T_span], h[t] == ϕ * w[t] * ΔT + B * z_on[t] * ΔT)

    # Bilinear term linearization (ω_t = A_t * z_t)
    @constraint(model, [t in T_span], w[t] <= α_max * z_on[t])
    @constraint(model, [t in T_span], w[t] >= α_min * z_on[t])
    @constraint(model, [t in T_span], w[t] >= A[t] - α_max * (1 - z_on[t]))
    @constraint(model, [t in T_span], w[t] <= A[t] - α_min * (1 - z_on[t]))

    # 3 state logic
    @constraint(model, [t in T_span], z_on[t] + z_off[t] == 1)

    # meet a daily demand quota
    n_hours_per_day = 24
    total_days = Int.(floor(length(T_span)/n_hours_per_day))
    D_H = 750 #sigma
    daily_indices = [
        T_span[(d-1)*n_hours_per_day + 1 : d*n_hours_per_day]
        for d in 1:total_days
    ]

    # 3. Build the constraint by iterating over your new groups
    @constraint(model, daily_quota[i in daily_indices],
        sum(h[t] for t in i) >= D_H
    )
    # # Start up binary logic
    @constraint(model, [t=2:length(T_span)], z_start[t] >= z_on[t] - z_on[t-1] ) # binary constraints
    @constraint(model, [t=2:length(T_span)], z_start[t] <= z_on[t])
    @constraint(model, [t=2:length(T_span)], z_start[t] <= 1 - z_on[t-1])

    # Economic-specific parameters
    λ_H = θ.λ_H # hydrogen selling price [$/kg]
    λ_CAPEX_yearly = θ.λ_CAPEX_yearly
    λ_CAPEX = θ.λ_CAPEX_Plant
    λ_OPEX = θ.λ_OPEX
    λ_CAPEX_replace = θ.λ_CAPEX_Stack

    @expression(model, cap_exp, λ_CAPEX * ϕ) # [=] $/MW * MW
    @expression(model, op_exp[y in years], λ_OPEX)
    @expression(model, replace_exp[y in years], ϕ * λ_CAPEX_replace * z_replace[y])
    # @expression(model, h_profit, sum(h[t] * λ_H for t in T_span)) # [=] kg H2 * $/kg H2
    @expression(model, h_revenue[y in year_indices], sum(h[t] * λ_H for t in y))
    @expression(model, e_exp[y in year_indices], sum(e_tot[t] * spp_avg_arr[t] for t in y)) # [=] MWhr * $/MWhr
    # @expression(model, start_costs, sum(z_start[t] * start_cost for t in T_span))
    # @expression(model, profit, h_profit - elec_exp - cap_exp - op_exp - replace_exp) # [=] $
    
    @expression(model, NPV[(y, y_k) in zip(years, year_indices)], (1/(1+ρ)^(y-1))* (-replace_exp[y] - e_exp[y_k] - op_exp[y] + h_revenue[y_k]))
    # Maximize profit
    @objective(model, Min, -sum(NPV[y] for y in zip(years, year_indices)))
    @expression(model, value, sum(NPV[y] for y in zip(years, year_indices)) - cap_exp)
    set_optimizer_attribute(model, "MIPGap", gap)
    # set_optimizer_attribute(model, "IterationLimit", 500000)
    # set_optimizer_attribute(model, "TimeLimit", 500)

    @expression(model, factor[y in years], (1/(1+ρ)^(y-1)))
    @expression(model, LCOH_num, λ_CAPEX + sum(factor[y]*(e_exp[y_k] + op_exp[y] + replace_exp[y]) for (y, y_k) in zip(years, year_indices)))
    @expression(model, LCOH_den, sum(factor[y]*sum(h[t] for t in y_k) for (y, y_k) in zip(years, year_indices)))
    @expression(model, LCOH, LCOH_num/LCOH_den)


    optimize!(model)

    return model
end