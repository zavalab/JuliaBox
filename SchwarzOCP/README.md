# SchwarzOCP

Overlapping Schwarz decomposition for nonlinear optimal control: convergence, effect of overlap, exponential decay, and quadrotor and thin-plate benchmarks.

## How to run
To reproduce the results, run the following script
```julia
julia --project=. -t 20 # make sure that your hardware has enough cores

] add https://github.com/sshin23/SimpleNLModels.jl#2a0023d802ba8211f60a8c920a1a5bbd18eb04ed
] add https://github.com/sshin23/SimpleSchwarz.jl#6ef555336f30ee2e6cedaf874e89f60a136cb2b2
] add https://github.com/sshin23/SimpleADMM.jl#408cb08d86fb76c8a08e07d0df9ac7d3241e1268
] instantiate

include("run_case_study.jl")
```

## Reference
This code accompanies the paper:

Sen Na, Sungho Shin, Mihai Anitescu, and Victor M. Zavala. "On the Convergence of Overlapping Schwarz Decomposition for Nonlinear Optimal Control." *IEEE Transactions on Automatic Control* 67, no. 11 (2022): 5996-6011. https://doi.org/10.1109/tac.2022.3194087

Preprint: https://arxiv.org/abs/2005.06674

## BibTeX
```bibtex
@article{na2022convergence,
  title   = {{On the Convergence of Overlapping Schwarz Decomposition for Nonlinear Optimal Control}},
  author  = {Na, Sen and Shin, Sungho and Anitescu, Mihai and Zavala, Victor M.},
  journal = {IEEE Transactions on Automatic Control},
  volume  = {67},
  number  = {11},
  pages   = {5996--6011},
  year    = {2022},
  doi     = {10.1109/tac.2022.3194087}
}
```

## Authors
Sen Na, Sungho Shin, Mihai Anitescu, and Victor M. Zavala
