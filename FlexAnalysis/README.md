# FlexAnalysis

Flexibility analysis case studies: a 2D example and 3-node, IEEE 14-bus, and 141-bus power networks.

## Contents
- `2d_example.jl`, `3node.jl`, `IEEE_14.jl`, `141node.jl`: case studies
- `data_reader.jl`, `141bus.txt`: data for the 141-bus network

These case studies use [FlexJuMP.jl](https://github.com/zavalab/FlexJuMP.jl).

## Reference
This code accompanies the paper:

Joshua L. Pulsipher, Daniel Rios, and Victor M. Zavala. "A computational framework for quantifying and analyzing system flexibility." *Computers & Chemical Engineering* 126 (2019): 342-355. https://doi.org/10.1016/j.compchemeng.2019.04.024

Preprint: https://arxiv.org/abs/2106.12706

## BibTeX
```bibtex
@article{pulsipher2019computational,
  title   = {{A computational framework for quantifying and analyzing system flexibility}},
  author  = {Pulsipher, Joshua L. and Rios, Daniel and Zavala, Victor M.},
  journal = {Computers \& Chemical Engineering},
  volume  = {126},
  pages   = {342--355},
  year    = {2019},
  doi     = {10.1016/j.compchemeng.2019.04.024}
}
```

## Authors
Joshua L. Pulsipher, Daniel Rios, and Victor M. Zavala
