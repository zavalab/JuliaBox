# AdchemCaseStudy

Graph-based modeling and decomposition of coupled natural gas and power networks using Plasmo.jl, with transient gas flow and multi-period optimal power flow case studies.

## Contents
- `AdchemCaseStudy.jl`: module that loads the models
- `ng_transient.jl`: transient natural gas network model (GasModels.jl)
- `mp_strg_opf.jl`: multi-period optimal power flow with storage (PowerModels.jl)
- `adchem_cstudy.jl`: runs the case studies; `adchem_plot.jl`: makes the figures
- `data/`: PGLib power network cases and gas network data and time series
- `output/`: figures

## Reference
This code accompanies the paper:

Sungho Shin, Carleton Coffrin, Kaarthik Sundar, and Victor M. Zavala. "Graph-Based Modeling and Decomposition of Energy Infrastructures." *IFAC-PapersOnLine* 54, no. 3 (2021): 693-698. https://doi.org/10.1016/j.ifacol.2021.08.322

Preprint: https://arxiv.org/abs/2010.02404

## BibTeX
```bibtex
@article{shin2021graph,
  title   = {{Graph-Based Modeling and Decomposition of Energy Infrastructures}},
  author  = {Shin, Sungho and Coffrin, Carleton and Sundar, Kaarthik and Zavala, Victor M.},
  journal = {IFAC-PapersOnLine},
  volume  = {54},
  number  = {3},
  pages   = {693--698},
  year    = {2021},
  doi     = {10.1016/j.ifacol.2021.08.322}
}
```

## Authors
Sungho Shin, Carleton Coffrin, Kaarthik Sundar, and Victor M. Zavala
