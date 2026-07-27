# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Julia research codebase that models the economics of hydrogen electrolyzers powered by grid
electricity, wind power, or both, under the IRA 45V clean hydrogen production tax credit. Each
configuration is formulated as a mixed-integer program (MIP) in JuMP, solved with Gurobi, over an
hourly time series (typically one year, 8760 hours).

## Running

There is no build/test/lint tooling — this is run interactively (e.g. in the Julia REPL or via
`julia --project=. main.jl`). Running requires:
- A licensed Gurobi installation (`Gurobi.Env()` is created at the top of `main.jl`).
- `Project.toml`/`Manifest.toml` are checked in; run `julia --project=. -e 'import Pkg; Pkg.instantiate()'`
  once to materialize the exact dependency set (`JuMP`, `Gurobi`, `DataFrames`, `XLSX`, `Plots`,
  `CPUTime`, `Statistics`, `Dates`, `JLD2`).
- The `data/` directory populated with the expected spreadsheets (see Data below). Only 2024 grid
  emissions data is currently checked in and required — 2022/2023/2025 are optional and are skipped
  with a `@warn` by `load_emissions()` if their source files aren't present under `data/<year>/`.

`main.jl` is the entry point: it loads emissions and day-ahead price data, builds a base-case
`Electrolyzer`, runs all three formulations (grid-only, wind-only, hybrid) back to back, and saves each
result to `results/basecase_{grid,wind,both}.jld2`. All paths inside `main.jl` and the `src/` files are
resolved relative to the file's own location (`@__DIR__`), so it can be launched from any working
directory.

## Architecture

Source lives in `src/`; `main.jl` at the repo root is the only script meant to be run directly.

### Core data type: `Electrolyzer` (`src/electrolyzer.jl`)

An immutable struct holding both raw device/economic parameters and several derived quantities computed
in its inner constructor (production coefficients `α_max`/`α_min` from efficiency and degradation
physics, capital recovery factor `CRF`, yearly/replacement CAPEX, etc.). It also carries a set of
`diff*` sensitivity-analysis multipliers (`diffcp`, `diffcr`, `diffeff`, `diffem`, `diffspp`, `difftier`)
that scale turbine efficiency, IRA credit value, electrolyzer efficiency, grid carbon intensity, power
price, and 45V emissions-tier thresholds respectively — these are the knobs used for scenario/sensitivity
sweeps rather than editing the model formulations. The constructor takes all 17 positional args
(including `difftier` last) — a previous version of the entry-point script omitted it, which is why it
no longer matches `main.jl`. This file is guarded with `if !@isdefined(Electrolyzer)` so it can be
safely `include`d multiple times from different scripts.

### Three optimization formulations, each with its own entry function

- `src/formulate_grid.jl` — `run_grid_opt(env, θ, spp_array, CI_grid; gap, z_fixed)`. Grid-only power
  source. Electrolyzer draws `e_grid`, incurs emissions `e_y = ci_grid .* e_grid`, pays the day-ahead
  settlement price `spp_array`.
- `src/formulate_wind.jl` — `run_wind_opt(env, θ, spp_array, CI_grid; gap, z_fixed)`. Wind-only power
  source. Includes `load_windspeeds.jl` to compute a 1 MW turbine's hourly available energy (`e_avail`,
  physics-based from wind speed) and models excess wind sold back to the grid (`e_sell`) at `spp_array`.
- `src/formulate_wind_grid.jl` — `run_hybrid_opt(env, θ, spp_array, CI_grid; gap, z_fixed)`. Hybrid:
  electrolyzer can draw from both `e_wind` and `e_grid` simultaneously (`e_tot = e_grid + e_wind`),
  combining the emissions accounting of the grid model with the wind-availability/sell-back logic of the
  wind model. This is the version `run_IRAcasestudy.jl` (now removed/superseded by `main.jl`) used to run.

These three formulations used to all define an identically-named `run_3st_opt`, which meant `include`ing
more than one in the same session silently made the last one win — a real bug, not just a naming
nit. They now have distinct names (`run_grid_opt`/`run_wind_opt`/`run_hybrid_opt`) specifically so
`main.jl` can run all three in one process. Keep that distinction if you add a fourth formulation.

All three share the same modeling pattern and are the natural place to look for how a given mechanic
(startup logic, tier credits, etc.) is implemented — check whether a fix belongs in one file or all three.

Common structure across all three formulations:
- **Three-state unit commitment**: each hour is exactly one of on / off / standby
  (`z_on + z_off + z_sb == 1`), with standby allowing a reduced power draw (`ρ_sb` fraction) and a
  `z_start` binary tracking cold starts (transitions from off/standby into on).
- **Hydrogen production**: `h[t] = α_max * (power drawn)[t] + B * z_on[t]`, where `B = 9.66` is a
  base/idle term and `α_max` is the (possibly `diffeff`-scaled) production coefficient from `θ`.
- **45V tiered tax credit**: implemented via a big-M formulation over 5 binary `tier` variables per hour
  rather than the commented-out `@constraint(... => {...})` indicator-constraint version left in the
  file for reference. Tier boundaries (`ϵ_upper`/`ϵ_lower`) and credit multipliers (`cr`) are scaled by
  `difftier`/`diffcr` from `θ`. Do not delete the commented indicator-constraint blocks without checking
  whether they document intended future behavior.
- **Objective**: maximize `Net_Profit` (hydrogen revenue + IRA credit + wind sell-back revenue − energy
  purchase cost − OPEX), implemented as `Min -Net_Profit`.
- Each function prints a solve summary (status, revenues, emissions, on/off/standby hour counts) inside
  a `try/catch` after `optimize!`, then returns the raw `model`.

### Data loading

- `src/load_emissions.jl` — reads ERCOT average-emissions-rate spreadsheets per year (2022–2025; 2022
  and 2023 use a different raw sheet layout than 2024/2025, hence separate parsing paths that both
  funnel into `clean_emissions`). `clean_emissions` builds a clean 8760-length hourly series per year,
  filling gaps via neighbor-hour averaging and handling the 2023 DST case (March 11 is missing and is
  interpolated as the average of March 10 and March 12). `load_emissions()` returns a
  `Dict{Int, Vector{Float64}}` keyed by year, only including years whose source files are present under
  `data/<year>/` — missing years are skipped with `@warn` rather than erroring. Returns emissions in the
  source units; callers typically rescale (e.g. `emissions_by_year[2024] .* 1000` in `main.jl`) to get
  `CI_grid`.
- `src/load_windspeeds.jl` — reads hourly 2024 wind speed data from `data/Data_weather.xlsx` (sheet
  `"2024"`, column 7), returns an 8760-length series consumed by the wind physics calc in the wind
  formulations.
- Day-ahead market prices (`spp_array` / `CI_grid` inputs) come from `data/ERCOT_DAM_AVG_2014-2024.xlsx`
  and `data/ERCOT_DAM_PAN_2014-2024.xlsx`, loaded directly in `main.jl` (15-minute interval data averaged
  down to hourly via `reshape`/`mean`).

### Persisting results (`src/save_results.jl`)

`populate_optimization_result` extracts an `OptimizationResult` (decision variable values, objective,
MIP gap, status, the `Electrolyzer` config, solve time, LCOH) from a solved JuMP `model` and
`pop_and_save` serializes it whole via `JLD2` (`jldopen`). `get_model_stats` reloads a saved result and
reports on/off/standby/cold-start/warm-start hour counts — use this rather than re-solving when just
inspecting a previous run. Saved results live in `results/` (`basecase_grid.jld2`, `basecase_wind.jld2`,
`basecase_both.jld2`), one per formulation, and are regenerated in place each time `main.jl` runs.

## Notes

- Greek-letter variable names (`ϕ`, `α`, `δ`, `η`, `ρ`, `θ`, `ξ`, `ℓ`, `λ`, `ϵ`) are used throughout to
  mirror the underlying engineering/economics notation; keep this convention when extending the models.
- No automated tests exist. Validate changes by running `main.jl` (or the relevant `run_*_opt` function
  directly) and comparing solver status/printed summary output — and, where relevant, `get_model_stats`
  on the saved result — against a known-good baseline in `results/`.
