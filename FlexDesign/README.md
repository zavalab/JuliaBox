# FlexDesign

Stochastic programming for the design of flexible systems, with 3-node, IEEE 14-bus, and 141-bus power network case studies.

## Contents
For each network (`3node`, `ieee14`, `141node`):
- `*_cont.jl`: design with continuous decisions
- `*_bin.jl`: design with binary decisions
- `*_analysis.jl`: flexibility analysis of the resulting designs
- `*_data*.json`: results

`data_reader.jl` and `141bus.txt` provide the data for the 141-bus network.

## Reference
This code accompanies the paper:

Joshua L. Pulsipher and Victor M. Zavala. "A scalable stochastic programming approach for the design of flexible systems." *Computers & Chemical Engineering* 128 (2019): 69-76. https://doi.org/10.1016/j.compchemeng.2019.05.033

Preprint: https://arxiv.org/abs/2106.12708

## BibTeX
```bibtex
@article{pulsipher2019scalable,
  title   = {{A scalable stochastic programming approach for the design of flexible systems}},
  author  = {Pulsipher, Joshua L. and Zavala, Victor M.},
  journal = {Computers \& Chemical Engineering},
  volume  = {128},
  pages   = {69--76},
  year    = {2019},
  doi     = {10.1016/j.compchemeng.2019.05.033}
}
```

## Authors
Joshua L. Pulsipher and Victor M. Zavala
