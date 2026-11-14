# WindSigVaR

Chance-constrained optimization of wind turbine pitch and torque control under wind-speed uncertainty using the sigmoidal value-at-risk (SigVaR) approximation, compared with CVaR, almost-sure, and smooth approximations.

## How to run
Run `task.jl` with Julia. It calls `AS.jl` (almost-sure), `CVar.jl` (CVaR), `svar.jl` (SigVaR), and `smooth.jl` (smooth approximation). Results are saved in the corresponding `.out` files (see `Result/`).

## Contents
- `setup.jl`, `createModel.jl`, `createModel_svar.jl`, `initial.jl`: wind turbine model, scenarios, and initialization
- `simulator.jl`, `plot.jl`: simulation and figures
- `ipopt.opt`, `pipsnlp.parameter`: solver options

## Reference
This code accompanies the paper:

Yankai Cao and Victor M. Zavala. "A Sigmoidal Approximation for Chance-Constrained Nonlinear Programs." arXiv preprint arXiv:2004.02402 (2020). https://doi.org/10.48550/arXiv.2004.02402

## BibTeX
```bibtex
@misc{cao2020sigmoidal,
  title         = {{A Sigmoidal Approximation for Chance-Constrained Nonlinear Programs}},
  author        = {Cao, Yankai and Zavala, Victor M.},
  year          = {2020},
  eprint        = {2004.02402},
  archivePrefix = {arXiv},
  doi           = {10.48550/arXiv.2004.02402}
}
```

## Authors
Yankai Cao and Victor M. Zavala
