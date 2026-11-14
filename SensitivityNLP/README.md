# SensitivityNLP

Exponential decay of sensitivity in graph-structured nonlinear programs: optimal control, optimal power flow, PDE-constrained, and stochastic examples with theoretical bounds.

## Contents
- `ocp.jl`, `opf.jl`, `pde.jl`, `stoch.jl`: optimal control, optimal power flow (PGLib 500-bus system), PDE-constrained, and stochastic examples
- `theoretical_bounds.jl`: theoretical bounds
- `data/`: data; `fig/`: figures
- `Project.toml`, `Manifest.toml`: Julia environment

## Reference
This code accompanies the paper:

Sungho Shin, Mihai Anitescu, and Victor M. Zavala. "Exponential Decay of Sensitivity in Graph-Structured Nonlinear Programs." *SIAM Journal on Optimization* 32, no. 2 (2022): 1156-1183. https://doi.org/10.1137/21m1391079

Preprint: https://arxiv.org/abs/2101.03067

## BibTeX
```bibtex
@article{shin2022exponential,
  title   = {{Exponential Decay of Sensitivity in Graph-Structured Nonlinear Programs}},
  author  = {Shin, Sungho and Anitescu, Mihai and Zavala, Victor M.},
  journal = {SIAM Journal on Optimization},
  volume  = {32},
  number  = {2},
  pages   = {1156--1183},
  year    = {2022},
  doi     = {10.1137/21m1391079}
}
```

## Authors
Sungho Shin, Mihai Anitescu, and Victor M. Zavala
