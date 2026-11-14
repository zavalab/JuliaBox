# PIPS_examples

Examples for PIPS-NLP: stochastic farmer problems and line-pack dispatch in a natural gas network, serial and parallel.

## Contents
- `farm_PIPS.jl`, `farmlarge_PIPS.jl`, `farm_PIPS*.ipynb`: stochastic farmer problems
- `stochpdegas_PIPS.jl`, `stochpdegas_PIPS*.ipynb`, `set_data_paper.jl`, `DAT_PAPER/`: stochastic natural gas network
- `ParPipsNlp.jl`, `NetParPipsNlp.jl`, `NetJuMP.jl`: interfaces to PIPS-NLP
- `pipsnlp.parameter`: PIPS-NLP options

## Reference
This code accompanies the paper:

Naiyuan Chiang, Cosmin G. Petra, and Victor M. Zavala. "Structured nonconvex optimization of large-scale energy systems using PIPS-NLP." *2014 Power Systems Computation Conference* (2014): 1-7. https://doi.org/10.1109/pscc.2014.7038374

## BibTeX
```bibtex
@article{chiang2014structured,
  title   = {{Structured nonconvex optimization of large-scale energy systems using PIPS-NLP}},
  author  = {Chiang, Naiyuan and Petra, Cosmin G. and Zavala, Victor M.},
  journal = {2014 Power Systems Computation Conference},
  pages   = {1--7},
  year    = {2014},
  doi     = {10.1109/pscc.2014.7038374}
}
```

## Authors
Naiyuan Chiang, Cosmin G. Petra, and Victor M. Zavala
