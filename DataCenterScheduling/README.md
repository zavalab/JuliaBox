# DataCenterScheduling

Receding-horizon scheduling model for data center operations that responds to power grid signals such as electricity prices and carbon emission rates.

## Contents
- `main.jl`: scheduling model; `run.jl`: runs the case study
- `data_gen.jl`: generates the workload data
- `0725_0801_data.csv`: grid data for July 25 to August 1

## Reference
This code accompanies the paper:

Weiqi Zhang, Line Roald, and Victor Zavala. "Exploring the Impacts of Power Grid Signals on Data Center Operations Using a Receding-Horizon Scheduling Model." *IEEE Transactions on Sustainable Computing* 8, no. 2 (2023): 245-256. https://doi.org/10.1109/tsusc.2022.3224668

Preprint: https://arxiv.org/abs/2204.06654

## BibTeX
```bibtex
@article{zhang2023exploring,
  title   = {{Exploring the Impacts of Power Grid Signals on Data Center Operations Using a Receding-Horizon Scheduling Model}},
  author  = {Zhang, Weiqi and Roald, Line and Zavala, Victor},
  journal = {IEEE Transactions on Sustainable Computing},
  volume  = {8},
  number  = {2},
  pages   = {245--256},
  year    = {2023},
  doi     = {10.1109/tsusc.2022.3224668}
}
```

## Authors
Weiqi Zhang, Line Roald, and Victor Zavala
