# Battery_FP

Stochastic model predictive control of a stationary battery in energy and frequency regulation markets while mitigating demand charges.

## Contents
- `setup.jl`: setup and data
- `FP_MPC.jl`, `MPC_flexible.jl`: MPC formulations
- `sim_MPC.jl`, `modified_sim_MPC.jl`: closed-loop simulations
- `Heuristic.jl`: heuristic policy for comparison

## Reference
This code accompanies the paper:

Ranjeet Kumar, Michael J. Wenzel, Matthew J. Ellis, Mohammad N. ElBsat, Kirk H. Drees, and Victor M. Zavala. "A Stochastic Model Predictive Control Framework for Stationary Battery Systems." *IEEE Transactions on Power Systems* 33, no. 4 (2018): 4397-4406. https://doi.org/10.1109/tpwrs.2017.2789118

## BibTeX
```bibtex
@article{kumar2018stochastic,
  title   = {{A Stochastic Model Predictive Control Framework for Stationary Battery Systems}},
  author  = {Kumar, Ranjeet and Wenzel, Michael J. and Ellis, Matthew J. and ElBsat, Mohammad N. and Drees, Kirk H. and Zavala, Victor M.},
  journal = {IEEE Transactions on Power Systems},
  volume  = {33},
  number  = {4},
  pages   = {4397--4406},
  year    = {2018},
  doi     = {10.1109/tpwrs.2017.2789118}
}
```

## Authors
Ranjeet Kumar, Michael J. Wenzel, Matthew J. Ellis, Mohammad N. ElBsat, Kirk H. Drees, and Victor M. Zavala
