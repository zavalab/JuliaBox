# FlareDesign

Flare system design under uncertainty using chance constraints, with conditional value-at-risk (CVaR) and sigmoidal value-at-risk (SVaR) approximations for normal, lognormal, exponential, and Gumbel distributions.

## Contents
- `*DistributionFlare.jl`, `MultiUncertaintyFlare.jl`: flare design models for each type of uncertainty
- `*DCvar.jl`, `*DCvarSvar.jl`: CVaR and CVaR+SVaR formulations (prefixes N, LN, E, and G stand for normal, lognormal, exponential, and Gumbel)
- `EDAS.jl`, `LogNormalAS.jl`: almost-sure formulations

## Reference
This code accompanies the paper:

Javier Tovar-Facio, Yankai Cao, Jose M. Ponce-Ortega, and Victor M. Zavala. "Scalable Solution Strategies for Chance-Constrained Nonlinear Programs." *Industrial & Engineering Chemistry Research* 57, no. 23 (2018): 7987-7998. https://doi.org/10.1021/acs.iecr.8b00646

## BibTeX
```bibtex
@article{tovarfacio2018scalable,
  title   = {{Scalable Solution Strategies for Chance-Constrained Nonlinear Programs}},
  author  = {Tovar-Facio, Javier and Cao, Yankai and Ponce-Ortega, Jose M. and Zavala, Victor M.},
  journal = {Industrial \& Engineering Chemistry Research},
  volume  = {57},
  number  = {23},
  pages   = {7987--7998},
  year    = {2018},
  doi     = {10.1021/acs.iecr.8b00646}
}
```

## Authors
Javier Tovar-Facio, Yankai Cao, Jose M. Ponce-Ortega, and Victor M. Zavala
