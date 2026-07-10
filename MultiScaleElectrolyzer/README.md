# Multi-Scale Optimization Framework for Grid-Integrated Electrolysis

Code to reproduce the results and figures from:

> K.X. Jennings, V.M. Zavala, S. Avraamidou. **"A Multi-Scale Optimization Framework for Grid-Integrated Electrolysis."** Preprint submitted to *Computers & Chemical Engineering*.

This repository implements a mixed-integer linear program (MILP) that co-optimizes short-term (hourly) demand response operation with long-term (annual) electrolyzer stack replacement decisions, accounting for usage-based degradation. The case study evaluates a 2.2 MW alkaline water electrolyzer participating in the ERCOT day-ahead market (DAM) over a 22-year horizon.

## Requirements

- **Julia** 1.12.5 (the `Manifest.toml` was generated with this version; other 1.1x versions will likely work but are untested)
- **Gurobi** ≥ 10.0 with a valid license (academic licenses are free — see [Gurobi's website](https://www.gurobi.com/academia/academic-program-and-licenses/)). The MILP formulation uses `Gurobi.jl` directly, so a Gurobi installation and license are **required** to run the optimization; other solvers are not currently supported.
- **Python 3** with `matplotlib`/`numpy` (only needed for the two `.py` plotting scripts in `plots/`)

Key Julia package versions (see `Manifest.toml` for the full dependency tree):

| Package | Version |
|---|---|
| JuMP | 1.30.0 |
| Gurobi | 1.9.2 |
| DataFrames | 1.8.1 |
| XLSX | 0.11.0 |
| Plots | 1.41.6 |
| CSV | 0.10.16 |
| JSON | 1.4.0 |
| CPUTime | 1.0.0 |

## Setup

This folder currently ships a `Manifest.toml` but no `Project.toml`. To install the exact dependency versions used to generate the paper's results:

```julia
julia --project=.
] instantiate
```

If `instantiate` complains about a missing `Project.toml`, add the packages listed above manually instead:

```julia
] add JuMP Gurobi DataFrames XLSX Plots CSV JSON CPUTime
```

Then make sure Gurobi.jl can find your Gurobi installation and license (see the [Gurobi.jl README](https://github.com/jump-dev/Gurobi.jl) if you haven't linked a local Gurobi install before).

## Repository structure

```
MultiScaleElectrolyzer/
├── data/
│   ├── ERCOT_DAM_AVG_2014-2024.xlsx   # ERCOT hub-average DAM prices, 2014–2024
│   └── ERCOT_DAM_PAN_2014-2024.xlsx   # ERCOT panhandle-node DAM prices, 2014–2024
├── src/
│   ├── electrolyzer_struct.jl         # Electrolyzer device/economic parameter struct
│   ├── formulation.jl                 # Core MILP formulation (JuMP + Gurobi)
│   ├── save_optimization.jl           # Helpers to extract/save solved model results
│   ├── run_casestudy.jl               # Case Study 1: ERCOT market participation
│   ├── run_sensitivity.jl             # Case Study 2: tornado sensitivity analysis
│   └── run_heatmap.jl                 # Case Study 2: NPV/LCOH parameter heatmap sweep
└── plots/
    ├── plot_price_series.jl           # Figure 3 (LMP price series + power consumption)
    ├── plot_efficiency_vs_price.jl    # Figure 5 (efficiency decay over time)
    ├── plot_replacements.ipynb        # Figure 6 (cumulative degradation)
    ├── electricity_markets.jl         # Market price summary statistics (Table 1)
    ├── case_study1.py                 # Figure 4 (stacked bar chart with LCOH)
    ├── plot_heatmap.py                # Figure 8 (NPV/LCOH heatmaps)
    └── figures/                       # Output figures (PDF)
```

## Running the case studies

All scripts are run from the `src/` directory (or with paths adjusted accordingly, since data paths are relative to each script's location via `@__DIR__`).

**Case Study 1 — ERCOT market participation** (Figures 3–6, Table 2):
```julia
julia src/run_casestudy.jl
```
This solves the full 22-year MILP for both the hub-average and panhandle price profiles and saves results as JSON to `data/test/`.

> **Known issue:** `run_casestudy.jl` currently references a file named `ERCOT_PAN_2014-2024.xlsx`, but the file in `data/` is named `ERCOT_DAM_PAN_2014-2024.xlsx`. Rename the file or fix the path in the script before running.

**Case Study 2 — Sensitivity analysis / tornado plot** (Figure 7, Tables 3–5):
```julia
julia src/run_sensitivity.jl
```
Several parameter perturbations (efficiency, on-degradation, startup degradation, stack cost, standby load) are defined but commented out in the script — uncomment the block(s) you need before running. Results save to `data/casestudy2/` or `data/tornado/`.

**Design tradeoff heatmap** (Figure 8):
```julia
julia src/run_heatmap.jl
```
Sweeps combinations of initial efficiency and operating degradation rate over an 11-year horizon. Edit the `efficiencies` and `op_degrads` arrays at the top of the script to control the sweep resolution. Results save to `data/heatmap/`.

Each run script solves to a specified MIP gap (default 1%, adjustable via the `gap` keyword to `run_3st_opt`/`run_2st_opt`) using Gurobi, and writes results to JSON via `pop_and_save` in `save_optimization.jl`.

## Generating figures

Once the case study JSON results exist, the plotting scripts in `plots/` regenerate the paper figures:

- `plot_price_series.jl`, `plot_efficiency_vs_price.jl`, `electricity_markets.jl` — run directly with `julia`
- `plot_replacements.ipynb` — open in Jupyter (with an IJulia kernel) or JupyterLab
- `case_study1.py`, `plot_heatmap.py` — run with `python3`, reading the JSON results produced above

Output figures are written to `plots/figures/`.

## Data

- `ERCOT_DAM_AVG_2014-2024.xlsx`: Hourly hub-average settlement point prices (SPP) from the ERCOT day-ahead market, 2014–2024.
- `ERCOT_DAM_PAN_2014-2024.xlsx`: Hourly panhandle-node SPP from the ERCOT day-ahead market, 2014–2024 (selected for proximity to wind generation).

Both series are replicated (concatenated with themselves) inside the run scripts to construct the 22-year study horizon.

## Citing this work

If you use this code, please cite the corresponding paper (citation details to be updated upon publication).

## Contact

For questions about the code, contact the corresponding author or open an issue on this repository.

---
*This README was generated with the help of Claude (Anthropic).*
