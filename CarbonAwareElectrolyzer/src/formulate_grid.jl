#grid-only
using JuMP, Gurobi
include("electrolyzer.jl")

"""
Hourly Core Optimization Model (DAM) - 1-Year Version with Soft Emissions Limits
Grid-only power source.
"""

function run_grid_opt(env::Gurobi.Env, θ::Electrolyzer, spp_array, CI_grid; gap=0.01, z_fixed=nothing)
    model = Model(() -> Gurobi.Optimizer(env))
    set_optimizer_attribute(model, "MIPGap", gap)

    # --- Structural Extraction & Parameters ---
    ϕ = θ.ϕ 
    α_max = θ.α_max
    α_init = α_max
    α_min = θ.α_min
    δ_on = θ.δ_on
    δ_start = θ.δ_start
    ρ = θ.i
    B = 9.66
    ρ_sb = θ.ρ_sb
    ΔT = 1

    n_per_hour = 4
    n_years = Int.(ceil(length(spp_array)/(8760*n_per_hour)))
    years=1:1
    ci_grid = Vector{Float64}(vec(collect(CI_grid)))


    diffem= θ.diffem
    diffeff= θ.diffeff
    diffspp= θ.diffspp
    diffcp= θ.diffcp
    diffcr= θ.diffcr
    difftier= θ.difftier

    α_max= diffeff .* α_max
    ci_grid = diffem .* ci_grid
    spp_avg_arr = diffspp .* spp_array
    T_span = collect(1:length(spp_avg_arr))

    # --- Variables ---
    @variables model begin
        e_y[T_span] >= 0
        e_grid[T_span] >= 0  
        h[T_span] >= 0
        z_sb[T_span], Bin 
        z_off[T_span], Bin 
        z_on[T_span], Bin
        z_start[T_span], Bin
        tier[T_span, 1:5], Bin
        IRA[T_span] >= 0 
    end

    #total emissions definition
    @constraint(model, [t in T_span], e_y[t] == ci_grid[t] * e_grid[t])
    
    # Device model constraints
    @constraint(model, [t in T_span], e_grid[t] <= ϕ * z_on[t] * ΔT + ϕ * z_sb[t] * ΔT * ρ_sb)
    @constraint(model, [t in T_span], e_grid[t] >= 0.2 * (ϕ * z_on[t] * ΔT) + ϕ * z_sb[t] * ΔT * ρ_sb)

    #capacity constraints
    @constraint(model, [t in T_span], h[t] == α_max * e_grid[t] + B * z_on[t] * ΔT)

    # 3 state logic
    @constraint(model, [t in T_span], z_on[t] + z_off[t] + z_sb[t] == 1)

    #start up
    @constraint(model, [t in 2:length(T_span)], z_start[t] >= z_on[t] - z_on[t-1] - z_sb[t-1]) 
    @constraint(model, [t in 2:length(T_span)], z_start[t] <= z_on[t])
    @constraint(model, [t in 2:length(T_span)], z_start[t] <= 1 - z_on[t-1] - z_sb[t-1])

    #big M-constraint method
    cr= diffcr .* [0, 0.6, 0.75, 1, 3]
    
    ϵ_upper= difftier .* [ 100, 4, 2.5, 1.5, 0.45]
    ϵ_lower= difftier .* [4, 2.5, 1.5, 0.45, 0]

    M=1e9
    
    @constraint(model, [t in T_span], sum(tier[t,i] for i in 1:5) == 1)
    @constraint(model, [t in T_span, i in 1:5], IRA[t] <= h[t] * cr[i] + M*(1-tier[t,i]))
    @constraint(model, [t in T_span, i in 1:5], IRA[t] >= h[t] * cr[i] - M*(1-tier[t,i]))

    @constraint(model, [t in T_span, i in 1:5], e_y[t] <= ϵ_upper[i] * h[t] + M*(1-tier[t,i]))
    @constraint(model, [t in T_span, i in 1:5], e_y[t] >= ϵ_lower[i] * h[t] - M*(1-tier[t,i]))


    #economics
    λ_H = θ.λ_H 
    λ_CAPEX = θ.λ_CAPEX_Plant
    λ_OPEX = θ.λ_OPEX

    @expression(model, h_revenue, sum(h[t] * λ_H for t in T_span))
    @expression(model, cap_exp, λ_CAPEX * ϕ) 
    @expression(model, op_exp, λ_OPEX)
    @expression(model, e_exp, sum(e_grid[t] * spp_array[t] for t in T_span)) 
    @expression(model, IRA_rev, sum(IRA[t] for t in T_span))

    @expression(model, Net_Profit, -e_exp - op_exp + h_revenue + IRA_rev)

    #run
    @objective(model, Min, -Net_Profit)
    optimize!(model)
    println("Solver Status Summary: ", termination_status(model))
try
        tot_h = sum(JuMP.value(h[t]) for t in T_span)
        tot_emissions = sum(JuMP.value(e_y[t]) for t in T_span)
        annual_epsilon = tot_h > 0 ? tot_emissions / tot_h : 0.0
        println("annual epsilon = ", annual_epsilon)
        println("H revenue = ", JuMP.value(model[:h_revenue]))
        println("energy cost = ", JuMP.value(model[:e_exp]))
        println("OPEX = ", JuMP.value(model[:op_exp]))

        #electricty usage
        e_grid_tot = sum(JuMP.value(model[:e_grid][t]) for t in T_span)
        tot_emission = sum(ci_grid[t] * JuMP.value(model[:e_grid][t]) for t in T_span)
        println("total grid power taken:", e_grid_tot)
        println("total emissions:", tot_emission)
        IRA_tot = sum(JuMP.value(model[:IRA][t]) for t in T_span)
        println("IRA", IRA_tot)

        # Use > 0.5 for binary variables
        hours_off = findall(t -> JuMP.value(model[:z_off][t]) > 0.5, T_span)
        println("Total hours off: ", length(hours_off))

        annual_operating_profit = JuMP.value(model[:Net_Profit])
        println("Annual Operating Profit: \$", annual_operating_profit)
        
    catch e
        println("Error during calculation: ", e)
    end
    return model
end
