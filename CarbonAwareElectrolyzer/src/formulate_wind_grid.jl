#hybrid wind and grid
using JuMP, Gurobi
include("electrolyzer.jl")
include("load_windspeeds.jl")

"""
Hourly Core Optimization Model (DAM) - 1-Year Version with Soft Emissions Limits
Hybrid model: electrolyzer can draw from wind and grid simultaneously.
"""

function run_hybrid_opt(env::Gurobi.Env, θ::Electrolyzer, spp_array, CI_grid; gap=0.01, z_fixed=nothing)
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
    ci_grid = Vector{Float64}(vec(collect(CI_grid)))

    years = 1:1

    diffem= θ.diffem
    diffeff= θ.diffeff
    diffspp= θ.diffspp
    diffcp= θ.diffcp
    diffcr= θ.diffcr
    difftier= θ.difftier

    α_max= diffeff .* α_max
    ci_grid = diffem .* ci_grid
    spp_avg_arr = diffspp .* spp_array
    #spp_avg_arr = spp_array
    T_span = collect(1:length(spp_avg_arr))

    # --- Variables ---
    @variables model begin
        e_y[T_span] >= 0
        e_wind[T_span] >= 0
        e_grid[T_span] >= 0  
        e_tot[T_span] >= 0
        h[T_span] >= 0
        z_sb[T_span], Bin 
        z_off[T_span], Bin 
        z_on[T_span], Bin
        z_start[T_span], Bin
        tier[T_span, 1:5], Bin
        IRA[T_span] >= 0 
        e_sell[T_span] >= 0
    end

    #wind turbine energy calc for 1 MW turbine
    ρ_e = 1.225  # kg/m3 density of air
    D = 61    # diameter of turbine
    A_e = 0.25 * pi * D^2  # area of single blade
    v_e = load_windspeeds()
    cp_e =  diffcp * 0.4   # efficiency
    #cp_e = 0.4
    e_avail = (0.5 * ρ_e * A_e * (v_e.^3) * cp_e) / 1e6

    #annual terms and total energy definition
    @constraint(model, [t in T_span], e_y[t] == ci_grid[t] * e_grid[t])
    @constraint(model, [t in T_span], e_tot[t] == e_grid[t] + e_wind[t])

    # Device model constraints (turned into flexible load problem)
    @constraint(model, [t in T_span], e_tot[t] <= ϕ * z_on[t] * ΔT + ϕ * z_sb[t] * ΔT * ρ_sb)
    @constraint(model, [t in T_span], e_tot[t] >= 0.2*(ϕ * z_on[t] * ΔT) + ϕ * z_sb[t] * ΔT * ρ_sb)

    #capacity constraints
    #@constraint(model, [t in T_span], h[t] == ϕ * w[t] * ΔT + B * z_on[t] * ΔT)
    @constraint(model, [t in T_span], e_wind[t] <= e_avail[t])  # Wind energy availability constraint
    @constraint(model, [t in T_span], e_sell[t] == e_avail[t] - e_wind[t])

    # 3 state logic
    @constraint(model, [t in T_span], z_on[t] + z_off[t] + z_sb[t] == 1)

    #start up
    @constraint(model, [t in 2:length(T_span)], z_start[t] >= z_on[t] - z_on[t-1] - z_sb[t-1]) 
    @constraint(model, [t in 2:length(T_span)], z_start[t] <= z_on[t])
    @constraint(model, [t in 2:length(T_span)], z_start[t] <= 1 - z_on[t-1] - z_sb[t-1])

    # hydrogen prod
    @constraint(model, [t in T_span], h[t] == e_tot[t] * ΔT * α_max + B * z_on[t] * ΔT)

    # Force  annual average emissions intensity to be clean 
    #@constraint(model, e_y[1] <= 4.0 * h_y[1])

    #indicator function method

    # Tier 1: No credit
    #@constraint(model, [y in years], tier[y, 1] => {e_y[y] >= 4.0 * h_y[y]})
    #@constraint(model, [y in years], tier[y, 1] => {IRA[y] == 0.0})

    # Tier 2: 0.60 credit
    #@constraint(model, [y in years], tier[y, 2] => {e_y[y] >= 2.5 * h_y[y]})
    #@constraint(model, [y in years], tier[y, 2] => {e_y[y] <= 4.0 * h_y[y]})
    #@constraint(model, [y in years], tier[y, 2] => {IRA[y] == 0.6 * h_y[y]})

    # Tier 3: 0.75 credit
    #@constraint(model, [y in years], tier[y, 3] => {e_y[y] >= 1.5 * h_y[y]})
    #@constraint(model, [y in years], tier[y, 3] => {e_y[y] <= 2.5 * h_y[y]})
    #@constraint(model, [y in years], tier[y, 3] => {IRA[y] == 0.75 * h_y[y]})

    # Tier 4: 1.00 credit
    #@constraint(model, [y in years], tier[y, 4] => {e_y[y] >= 0.45 * h_y[y]})
    #@constraint(model, [y in years], tier[y, 4] => {e_y[y] <= 1.5 * h_y[y]})
    #@constraint(model, [y in years], tier[y, 4] => {IRA[y] == 1.0 * h_y[y]})

    # Tier 5: 3.00 credit
    #@constraint(model, [y in years], tier[y, 5] => {e_y[y] <= 0.45 * h_y[y]})
    #@constraint(model, [y in years], tier[y, 5] => {IRA[y] == 3.0 * h_y[y]})

    #big M-constraint method
    #diffcr = θ.diffcr
    cr= diffcr .* [0, 0.6, 0.75, 1, 3]
    #cr = [0, 0.6, 0.75, 1, 3]
    
    
    ϵ_upper= difftier .* [ 100, 4, 2.5, 1.5, 0.45]
    ϵ_lower= difftier .* [4, 2.5, 1.5, 0.45, 0]

    M=1e9
    @constraint(model, [t in T_span], sum(tier[t,i] for i in 1:5) == 1)
    @constraint(model, [t in T_span, i in 1:5], IRA[t] <= h[t] * cr[i] + M*(1-tier[t,i]))
    @constraint(model, [t in T_span, i in 1:5], IRA[t] >= h[t] * cr[i] - M*(1-tier[t,i]))

    @constraint(model, [t in T_span, i in 1:5], e_y[t] <= ϵ_upper[i] * h[t] + M*(1-tier[t,i]))
    @constraint(model, [t in T_span, i in 1:5], e_y[t] >= ϵ_lower[i] * h[t] - M*(1-tier[t,i]))


    #economic
    λ_H = θ.λ_H 
    λ_CAPEX = θ.λ_CAPEX_Plant
    λ_OPEX = θ.λ_OPEX
   # level_cost = 42.0   # USD/MWh

    @expression(model, h_revenue, sum(h[t] * λ_H for t in T_span))
    @expression(model, cap_exp, λ_CAPEX * ϕ) 
    @expression(model, op_exp, λ_OPEX)
    @expression(model, e_exp, sum(e_grid[t] * spp_avg_arr[t] for t in T_span)) 
    #@expression(model, e_exp2, sum(e_wind[t] * level_cost for t in T_span))
    @expression(model, e_rev, sum(e_sell[t] * spp_avg_arr[t] for t in T_span))
    @expression(model, IRA_rev, sum(IRA[t] for t in T_span))

    @expression(model, Net_Profit, -e_exp - op_exp  + e_rev + h_revenue + IRA_rev)

    #run
    @objective(model, Min, -Net_Profit)
    optimize!(model)

    println("Total Variables: ", num_variables(model))
    println("Binary Variables: ", count(is_binary, all_variables(model)))
    # Replace complex calls using `count_variable_in_set_constraints` with:
    total_constraints = sum(
    num_constraints(model, F, S) 
    for (F, S) in list_of_constraint_types(model)
    )
    println("Total Constraints: ", total_constraints)
    println("Solver Status Summary: ", termination_status(model))
try
        tot_h = sum(JuMP.value(h[t]) for t in T_span)
        tot_emissions = sum(JuMP.value(e_y[t]) for t in T_span)
        annual_epsilon = tot_h > 0 ? tot_emissions / tot_h : 0.0
        println("annual epsilon = ", annual_epsilon)
        IRA_tot = sum(JuMP.value(model[:IRA][t]) for t in T_span)
        println("IRA", IRA_tot)
        println("H revenue = ", JuMP.value(model[:h_revenue]))
        println("energy cost = ", JuMP.value(model[:e_exp]))
        println("energy profit =", JuMP.value(model[:e_rev]))
        println("OPEX = ", JuMP.value(model[:op_exp]))

        #electricty usage
        e_wind_tot = sum(JuMP.value(model[:e_wind][t]) for t in T_span)
        e_grid_tot = sum(JuMP.value(model[:e_grid][t]) for t in T_span)
        tot_emission = sum(ci_grid[t] * JuMP.value(model[:e_grid][t]) for t in T_span)
        println("total wind power taken:", e_wind_tot)
        println("total grid power taken:", e_grid_tot)
        println("total emissions:", tot_emission)

        # Use > 0.5 for binary variables
        hours_off = findall(t -> JuMP.value(model[:z_off][t]) > 0.5, T_span)
        hours_on = findall(t -> JuMP.value(model[:z_on][t]) > 0.5, T_span)
        println("Total hours off: ", length(hours_off))
        println("Total hours on: ", length(hours_on))

        #println("max possible wind H2 left ≈ ", (sum(e_avail) - sum(JuMP.value(model[:e_wind][t]) for t in T_span)) * α_max)

        #println("H2 needed for tier 4 = ", JuMP.value(model[:e_y][1]) / 4.0 - JuMP.value(model[:h_y][1]))

        annual_operating_profit = JuMP.value(model[:Net_Profit])
        println("Annual Operating Profit: \$", annual_operating_profit)

        #can calc NPV
        #total_capex = θ.λ_CAPEX_Plant * θ.ϕ  
        #lifetime_years = 10
        #discount_rate = 0.05

        #present_value_of_revenues = sum(annual_operating_profit / (1 + discount_rate)^t for t in 1:lifetime_years)

        #project_NPV = -total_capex + present_value_of_revenues
        #println("True Project NPV: \$", project_NPV)
    catch e
        println("Error during calculation: ", e)
    end
    return model
end