# StochMOO

Scalable modeling and solution of stochastic multiobjective optimization problems: a combined heat and power (CHP) system with cost, emissions, and water objectives over 188 scenarios, solved with Plasmo and Ipopt.

## Contents
- `FormulationC.jl`, `FormulationD.jl`: stochastic multiobjective formulations
- `Utopian.jl`, `chp_model_up.jl`: utopia points
- `chp_model.jl`, `data.jl`: CHP scenario model and data
- `plotData.jl`, `MatlabFigures/`: results and figures

## Reference
This code accompanies the paper:

Yankai Cao, Luis Fabian Fuentes-Cortes, Siyu Chen, and Victor M. Zavala. "Scalable modeling and solution of stochastic multiobjective optimization problems." *Computers & Chemical Engineering* 99 (2017): 185-197. https://doi.org/10.1016/j.compchemeng.2017.01.021

## BibTeX
```bibtex
@article{cao2017scalable,
  title   = {{Scalable modeling and solution of stochastic multiobjective optimization problems}},
  author  = {Cao, Yankai and Fuentes-Cortes, Luis Fabian and Chen, Siyu and Zavala, Victor M.},
  journal = {Computers \& Chemical Engineering},
  volume  = {99},
  pages   = {185--197},
  year    = {2017},
  doi     = {10.1016/j.compchemeng.2017.01.021}
}
```

## Authors
Yankai Cao, Luis Fabian Fuentes-Cortes, Siyu Chen, and Victor M. Zavala
