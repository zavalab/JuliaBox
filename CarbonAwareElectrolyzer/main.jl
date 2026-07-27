using JuMP, Gurobi, DataFrames, XLSX, CPUTime, Statistics

# Resolve all paths relative to this file's location so `main.jl` can be run
# from any working directory (e.g. `julia --project=. main.jl` from anywhere).
const ROOT_DIR = @__DIR__
const DATA_DIR = joinpath(ROOT_DIR, "data")
const RESULTS_DIR = joinpath(ROOT_DIR, "results")
mkpath(RESULTS_DIR)

include(joinpath(ROOT_DIR, "src", "electrolyzer.jl"))
include(joinpath(ROOT_DIR, "src", "load_emissions.jl"))
include(joinpath(ROOT_DIR, "src", "load_windspeeds.jl"))
include(joinpath(ROOT_DIR, "src", "formulate_grid.jl"))
include(joinpath(ROOT_DIR, "src", "formulate_wind.jl"))
include(joinpath(ROOT_DIR, "src", "formulate_wind_grid.jl"))
include(joinpath(ROOT_DIR, "src", "save_results.jl"))

const env = Gurobi.Env()

# --- 1. Load input data ---
emissions_by_year = load_emissions()
haskey(emissions_by_year, 2024) || error("2024 grid emissions data is required (data/2024/) but was not found.")
CI_grid = emissions_by_year[2024] .* 1000

DAM_10yr = DataFrame(XLSX.readtable(joinpath(DATA_DIR, "ERCOT_DAM_AVG_2014-2024.xlsx"), "Sheet1"))
spp_arr_DAM = Float64.(DAM_10yr[:, "SettlementPointPrice"])
per_year = 35040
DAM_2024 = Float64.(spp_arr_DAM[1:per_year])

# 15-minute interval -> hourly
DAM_hr_2024 = vec(mean(reshape(DAM_2024, 8760, 4), dims=2))

# --- 2. Base-case electrolyzer parameters ---
θ_basecase = Electrolyzer(
    2.2,          # MW Capacity
    0.65,         # % LHV Efficiency
    3.2/1000000,  # V/hr degradation
    1.9,          # Minimum Volts
    10,           # Lifetime (years)
    0.0204/500,   # Start-up degradation
    3,            # $/kg H2 price
    0.05,         # Discount Rate
    1816*1000,    # Plant CAPEX
    250*1000,     # Stack CAPEX
    0.05,         # Standby load fraction
    1.0,          # efficiency of turbine
    1.0,          # credit multiplier
    1.0,          # max efficiency of electrolyzer
    1.0,          # ci_grid multiplier
    1.0,          # spp_array multiplier
    1.0,          # tier level diff scalar
)

# --- 3. Run all three case studies and save results ---
println("--- LAUNCHING OPTIMIZATION WORKSPACE ---")

println("\n[1/3] Grid-only case study")
t_grid = @CPUelapsed model_grid = run_grid_opt(env, θ_basecase, DAM_hr_2024, CI_grid)
pop_and_save(model_grid, t_grid, θ_basecase, joinpath(RESULTS_DIR, "basecase_grid.jld2"))

println("\n[2/3] Wind-only case study")
t_wind = @CPUelapsed model_wind = run_wind_opt(env, θ_basecase, DAM_hr_2024, CI_grid)
pop_and_save(model_wind, t_wind, θ_basecase, joinpath(RESULTS_DIR, "basecase_wind.jld2"))

println("\n[3/3] Hybrid wind + grid case study")
t_both = @CPUelapsed model_both = run_hybrid_opt(env, θ_basecase, DAM_hr_2024, CI_grid)
pop_and_save(model_both, t_both, θ_basecase, joinpath(RESULTS_DIR, "basecase_both.jld2"))

println("\n--- DONE. Results saved to $RESULTS_DIR ---")
