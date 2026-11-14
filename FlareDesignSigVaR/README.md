# FlareDesignSigVaR

Flare design with chance constraints using the sigmoidal value-at-risk (SigVaR) approximation, compared with almost-sure and smooth approximations.

## Contents
- `EDsvar.jl`: SigVaR approximation (the first iteration uses a CVaR approximation)
- `EDAS.jl`: almost-sure approach
- `EDsmooth.jl`: smooth approximation
- `*.out`: results; `ipopt.opt`: Ipopt options

See also `README`.

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
