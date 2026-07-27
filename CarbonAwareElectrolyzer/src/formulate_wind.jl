#wind-only
using JuMP, Gurobi
include("electrolyzer.jl")
include("load_windspeeds.jl")

"""
Hourly Core Optimization Model (DAM) - 1-Year Version with Soft Emissions Limits
Wind-only power source (1 MW turbine physics model), with excess wind sold to the grid.
"""

function run_wind_opt(env::Gurobi.Env, θ::Electrolyzer, spp_array, CI_grid; gap=0.01, z_fixed=nothing)
    model = Model(() -> Gurobi.Optimizer(env))
    set_optimizer_attribute(model, "MIPGap", gap)
    
    # Structural Extraction
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
    CRF = θ.CRF
    i = θ.i

    n_per_hour = 4
    n_years = Int.(ceil(length(spp_array)/(8760*n_per_hour)))
    ci_grid = Vector{Float64}(vec(collect(CI_grid)))
    years=1:1

    diffem= θ.diffem
    diffeff= θ.diffeff
    diffspp= θ.diffspp
    diffcp= θ.diffcp
    diffcr= θ.diffcr
    difftier= θ.difftier

    lifetime = 20
    cap_turbine = 1.125e6   #usd
    
    #wind turbine energy calc for 1 MW turbine
    ρ_e = 1.225  # kg/m3 density of air
    D = 61    # diameter of turbine
    A_e = 0.25 * pi * D^2  # area of single blade
    v_e = load_windspeeds()
    cp_e =  diffcp * 0.4   # efficiency
    #cp_e = 0.4
    e_avail = (0.5 * ρ_e * A_e * (v_e.^3) * cp_e) / 1e6

    α_max= diffeff .* α_max
    ci_grid = diffem .* ci_grid
    spp_avg_arr = diffspp .* spp_array
    #spp_avg_arr = spp_array
    T_span = collect(1:length(spp_avg_arr))

    # --- Variables ---
    @variables model begin 
        e_wind[T_span] >= 0
        h_tot[T_span] >= 0
        h[T_span] >= 0
        z_sb[T_span], Bin 
        z_off[T_span], Bin 
        z_start[T_span], Bin 
        z_on[T_span], Bin 
        e_sell[T_span] >= 0
    end


    # Operational Constraints
    @constraint(model, [t in T_span], e_wind[t] <= ϕ * z_on[t] * ΔT + ϕ * z_sb[t] * ΔT * ρ_sb)
    @constraint(model, [t in T_span], e_wind[t] >= 0.2*(ϕ * z_on[t] * ΔT) + ϕ * z_sb[t] * ΔT * ρ_sb)
    @constraint(model, [t in T_span], z_on[t] + z_off[t] + z_sb[t] == 1)

    #capacity constraints 
    @constraint(model, [t in T_span], e_wind[t] <= e_avail[t] * z_on[t])  # Wind energy availability constraint
    @constraint(model, [t in T_span], e_sell[t] == e_avail[t] - e_wind[t])

    #start up
    @constraint(model, [t in 2:length(T_span)], z_start[t] >= z_on[t] - z_on[t-1] - z_sb[t-1]) 
    @constraint(model, [t in 2:length(T_span)], z_start[t] <= z_on[t])
    @constraint(model, [t in 2:length(T_span)], z_start[t] <= 1 - z_on[t-1] - z_sb[t-1])

    #h prod
    @constraint(model, [t in T_span], h[t] == e_wind[t] * α_max + B * z_on[t] * ΔT)

    # 
    λ_H = θ.λ_H 
    λ_CAPEX = θ.λ_CAPEX_Plant
    λ_OPEX = θ.λ_OPEX
    level_cost=42   #USD/MWH
   
    cr= diffcr *3

    @expression(model, h_revenue, sum(h[t] * λ_H for t in T_span))
    @expression(model, cap_exp, λ_CAPEX * ϕ) 
    @expression(model, op_exp, λ_OPEX)
    #@expression(model, e_exp, sum(e_wind[t] * level_cost for t in T_span)) 
    @expression(model, IRA, sum(h[t] * cr for t in T_span))
    @expression(model, e_rev, sum(e_sell[t] * spp_array[t] for t in T_span))
    

    @expression(model, Net_Profit,  - op_exp + e_rev + h_revenue + IRA)
    
    # Objective 
    @objective(model, Min, -Net_Profit)

    # LCOH 
    @expression(model, LCOH, ((λ_CAPEX * ϕ)  + op_exp) / sum(h_tot[t] for t in T_span))

    optimize!(model)
    println("Solver Status Summary: ", termination_status(model))

try
        wind_prod = sum(JuMP.value(model[:h][t]) for t in T_span)
        println("Total wind farm production:", wind_prod)
        hours_wind = findall(t -> JuMP.value(model[:z_on][t]) > 0.5, T_span)
        println("Hours wind farm used ($(length(hours_wind)) hours):")
        hours_off = findall(t -> JuMP.value(model[:z_off][t]) > 0.5, T_span)
        println("Hours off ($(length(hours_off)) hours):")
        e_wind_tot = sum(JuMP.value(model[:e_wind][t]) for t in T_span)
        println("total wind power taken:", e_wind_tot)
        e_revenue= sum(JuMP.value(model[:e_sell][t]) * spp_array[t] for t in T_span)
        println("total revenue from selling electricity to grid:", e_revenue)

        println("H revenue = ", JuMP.value(model[:h_revenue]))
        println("OPEX = ", JuMP.value(model[:op_exp]))
        println("OPEX = ", JuMP.value(model[:cap_exp]))
        println("IRA = ", JuMP.value(model[:IRA]))

        annual_operating_profit = JuMP.value(model[:Net_Profit])
        println("Annual Operating Profit: \$", annual_operating_profit)

        #NPV calculation
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