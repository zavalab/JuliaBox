# CI Electrolyzer

Mixed-integer optimization models for the techno-economics of a PEM electrolyzer operating under the
IRA Section 45V clean hydrogen production tax credit, powered by grid electricity, wind, or both. Models
are formulated hourly (one year, 8760 hours) in [JuMP](https://jump.dev/) and solved with
[Gurobi](https://www.gurobi.com/).

## Requirements

- [Julia](https://julialang.org/) 1.10+ (developed against 1.12).
- A licensed Gurobi installation, reachable by `Gurobi.jl` (an academic license is sufficient).
- The input spreadsheets under `data/` (see [Data](#data) below).

## Setup

From the repository root:

```bash
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

This reads `Project.toml`/`Manifest.toml` and installs the exact dependency versions used to develop
the models (`JuMP`, `Gurobi`, `DataFrames`, `XLSX`, `Plots`, `CPUTime`, `JLD2`, `Statistics`, `Dates`).

## Running

```bash
julia --project=. main.jl
```

`main.jl` loads the emissions and day-ahead price data, builds a base-case `Electrolyzer` configuration,
and runs all three formulations back to back:

1. **Grid-only** (`src/formulate_grid.jl`, `run_grid_opt`)
2. **Wind-only** (`src/formulate_wind.jl`, `run_wind_opt`)
3. **Hybrid wind + grid** (`src/formulate_wind_grid.jl`, `run_hybrid_opt`)

Each run prints a solver summary (status, revenues, emissions, on/off/standby hour counts) and saves a
serialized `OptimizationResult` to `results/basecase_grid.jld2`, `results/basecase_wind.jld2`, and
`results/basecase_both.jld2` respectively, via `src/save_results.jl`. All file paths are resolved
relative to the repository root, so `main.jl` can be launched from any working directory.

To inspect a saved result without re-solving:

```julia
include("src/save_results.jl")
get_model_stats("results/basecase_both.jld2")
```

To run a single formulation with different parameters (e.g. a sensitivity sweep over the `diff*`
multipliers on `Electrolyzer`), include the relevant `src/formulate_*.jl` file and call its `run_*_opt`
function directly instead of going through `main.jl`.

## Data

Input spreadsheets live under `data/`:

- `data/<year>/` — ERCOT average-emissions-rate (AER) monthly reports, used by
  `src/load_emissions.jl`/`load_emissions()`. Only `data/2024/` is required — `main.jl` uses 2024 grid
  emissions as `CI_grid`. `2022`, `2023`, and `2025` are optional; if their source files aren't present,
  `load_emissions()` skips that year with a `@warn` instead of erroring.
- `data/Data_weather.xlsx` — hourly 2024 wind speed data, used by `src/load_windspeeds.jl`.
- `data/ERCOT_DAM_AVG_2014-2024.xlsx`, `data/ERCOT_DAM_PAN_2014-2024.xlsx` — ERCOT day-ahead market
  settlement point prices, loaded directly in `main.jl`.

## Repository structure

```
main.jl                     entry point — runs all three case studies
Project.toml, Manifest.toml Julia environment
src/
  electrolyzer.jl            Electrolyzer struct (device + economic parameters)
  formulate_grid.jl           grid-only MIP formulation
  formulate_wind.jl           wind-only MIP formulation
  formulate_wind_grid.jl      hybrid wind + grid MIP formulation
  load_emissions.jl           ERCOT emissions data loading/cleaning
  load_windspeeds.jl          wind speed data loading
  save_results.jl             JLD2 result serialization/inspection
data/                        input spreadsheets
results/                     saved OptimizationResult .jld2 files
```

## Mathematical formulation: hybrid wind + grid model

The formulation below is `run_hybrid_opt` in `src/formulate_wind_grid.jl`, the model used by `main.jl`.
The grid-only and wind-only formulations are restrictions of this model (drop the wind or grid terms
respectively).

### Sets

- $t \in \mathcal{T} = \{1, \dots, T\}$: hourly time steps ($T = 8760$ for one year).
- $i \in \{1, \dots, 5\}$: 45V emissions-intensity credit tiers.

### Parameters

| Symbol | Description |
|---|---|
| $\phi$ | Electrolyzer nominal capacity (MW) |
| $\alpha_{\max}$ | Hydrogen production coefficient (kg H$_2$/MWh) |
| $B = 9.66$ | Base/idle hydrogen production term while on |
| $\rho_{sb}$ | Standby power draw, as a fraction of $\phi$ |
| $\Delta T = 1$ | Time step length (hours) |
| $c_t$ | Grid carbon intensity at hour $t$ (kg CO$_2$/MWh) |
| $\pi_t$ | Day-ahead settlement point price at hour $t$ (\$/MWh) |
| $\bar{e}_t$ | Available wind energy at hour $t$ (MWh), from turbine physics |
| $\lambda_H$ | Hydrogen price (\$/kg) |
| $\lambda_{CAPEX}$ | Plant CAPEX rate (\$/MW) |
| $\lambda_{OPEX}$ | Annual O&M cost (\$) |
| $cr_i$ | 45V credit value for tier $i$ (\$/kg H$_2$): $[0,\ 0.60,\ 0.75,\ 1.00,\ 3.00]$ |
| $\epsilon^u_i,\ \epsilon^l_i$ | Upper/lower emissions-intensity bounds (kg CO$_2$/kg H$_2$) defining tier $i$ |
| $M$ | Big-M constant ($10^9$) |

$\bar{e}_t$ is computed from a 1 MW reference turbine: with air density $\rho_e = 1.225\ \text{kg/m}^3$,
blade diameter $D = 61\ \text{m}$, swept area $A_e = \tfrac{\pi}{4}D^2$, hub-height wind speed $v_t$, and
power coefficient $c_p$,

$$
\bar{e}_t = \frac{0.5\,\rho_e\,A_e\,v_t^3\,c_p}{10^6} .
$$

All parameters that vary in the code's sensitivity sweeps ($\alpha_{\max}$, $c_t$, $\pi_t$, $c_p$,
$cr_i$, $\epsilon^u_i, \epsilon^l_i$) are the base values scaled by an `Electrolyzer`-level multiplier
(`diffeff`, `diffem`, `diffspp`, `diffcp`, `diffcr`, `difftier`).

### Decision variables

| Symbol | Domain | Description |
|---|---|---|
| $e^{grid}_t$ | $\geq 0$ | Grid electricity drawn (MWh) |
| $e^{wind}_t$ | $\geq 0$ | Wind electricity drawn by the electrolyzer (MWh) |
| $e^{tot}_t$ | $\geq 0$ | Total electricity drawn (MWh) |
| $e^{sell}_t$ | $\geq 0$ | Excess wind energy sold to the grid (MWh) |
| $e^{y}_t$ | $\geq 0$ | Emissions attributed to hour $t$ (kg CO$_2$) |
| $h_t$ | $\geq 0$ | Hydrogen produced (kg) |
| $z^{on}_t, z^{off}_t, z^{sb}_t$ | $\{0,1\}$ | On / off / standby state indicators |
| $z^{start}_t$ | $\{0,1\}$ | Cold-start indicator |
| $\text{tier}_{t,i}$ | $\{0,1\}$ | 1 if hour $t$ falls in 45V tier $i$ |
| $IRA_t$ | $\geq 0$ | 45V credit revenue accrued at hour $t$ (\$) |

### Objective

Maximize annual net profit (implemented as $\min -\text{Net Profit}$):

$$
\max \; \text{Net Profit} = \underbrace{\sum_{t} \lambda_H\, h_t}_{h_{\text{revenue}}}
\;+\; \underbrace{\sum_{t} IRA_t}_{IRA_{\text{rev}}}
\;+\; \underbrace{\sum_{t} \pi_t\, e^{sell}_t}_{e_{\text{rev}}}
\;-\; \underbrace{\sum_{t} \pi_t\, e^{grid}_t}_{e_{\text{exp}}}
\;-\; \lambda_{OPEX}
$$

(Plant CAPEX, $\lambda_{CAPEX}\,\phi$, is computed as `cap_exp` but is not currently subtracted from the
objective — it is reported alongside the solve, not optimized against.)

### Constraints

**Emissions and energy balance**
$$
e^{y}_t = c_t \, e^{grid}_t \qquad \forall t
$$
$$
e^{tot}_t = e^{grid}_t + e^{wind}_t \qquad \forall t
$$

**Electrolyzer capacity (on/standby power envelope, minimum 20% turndown while on)**
$$
0.2\,\phi\, z^{on}_t\, \Delta T + \phi\, z^{sb}_t\, \Delta T\, \rho_{sb}
\;\leq\; e^{tot}_t \;\leq\;
\phi\, z^{on}_t\, \Delta T + \phi\, z^{sb}_t\, \Delta T\, \rho_{sb} \qquad \forall t
$$

**Wind availability and sell-back**
$$
e^{wind}_t \leq \bar{e}_t, \qquad e^{sell}_t = \bar{e}_t - e^{wind}_t \qquad \forall t
$$

**Three-state unit commitment**
$$
z^{on}_t + z^{off}_t + z^{sb}_t = 1 \qquad \forall t
$$

**Cold-start logic** ($t \geq 2$; $z^{start}_t = 1$ exactly when the unit turns on from off/standby)
$$
z^{start}_t \geq z^{on}_t - z^{on}_{t-1} - z^{sb}_{t-1}
$$
$$
z^{start}_t \leq z^{on}_t
$$
$$
z^{start}_t \leq 1 - z^{on}_{t-1} - z^{sb}_{t-1}
$$

**Hydrogen production**
$$
h_t = \alpha_{\max}\, e^{tot}_t\, \Delta T + B\, z^{on}_t\, \Delta T \qquad \forall t
$$

**45V tier assignment** (exactly one tier per hour)
$$
\sum_{i=1}^{5} \text{tier}_{t,i} = 1 \qquad \forall t
$$

**Tier-conditioned emissions-intensity bounds** (big-M: only active when $\text{tier}_{t,i}=1$)
$$
\epsilon^l_i\, h_t - M(1-\text{tier}_{t,i}) \;\leq\; e^{y}_t \;\leq\; \epsilon^u_i\, h_t + M(1-\text{tier}_{t,i})
\qquad \forall t,\ \forall i
$$

**Tier-conditioned credit value** (big-M)
$$
cr_i\, h_t - M(1-\text{tier}_{t,i}) \;\leq\; IRA_t \;\leq\; cr_i\, h_t + M(1-\text{tier}_{t,i})
\qquad \forall t,\ \forall i
$$

### 45V tiers

| Tier $i$ | Emissions intensity $\epsilon$ (kg CO$_2$/kg H$_2$) | Credit $cr_i$ (\$/kg H$_2$) |
|---|---|---|
| 1 | $\epsilon \geq 4.0$ | \$0.00 |
| 2 | $2.5 \leq \epsilon < 4.0$ | \$0.60 |
| 3 | $1.5 \leq \epsilon < 2.5$ | \$0.75 |
| 4 | $0.45 \leq \epsilon < 1.5$ | \$1.00 |
| 5 | $\epsilon < 0.45$ | \$3.00 |

The grid-only formulation (`src/formulate_grid.jl`) is this model with $e^{wind}_t = e^{sell}_t = 0$ and
$e^{tot}_t = e^{grid}_t$. The wind-only formulation (`src/formulate_wind.jl`) is this model with
$e^{grid}_t = 0$, $e^{tot}_t = e^{wind}_t$, and $e^y_t \equiv 0$ (no emissions accounting, since it never
draws from the grid).
