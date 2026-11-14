# PlasmoExamples

Example scripts and case studies for the OptiGraph abstraction in Plasmo.jl, including graph partitioning and decomposition of a natural gas network and a DC power grid.

## Contents
- `example1.jl` to `example5.jl`, `examples/`, `run_examples.jl`: introductory OptiGraph examples
- `case_studies/gasnetwork/`, `run_casestudy1.jl`: natural gas network case study
- `case_studies/dcpowergrid/`, `run_casestudy2.jl`: DC optimal power flow case study (PGLib 9241-bus system)
- `Project.toml`, `Manifest.toml`: Julia environment

## Reference
This code accompanies the paper:

Jordan Jalving, Sungho Shin, and Victor M. Zavala. "A graph-based modeling abstraction for optimization: concepts and implementation in Plasmo.jl." *Mathematical Programming Computation* 14, no. 4 (2022): 699-747. https://doi.org/10.1007/s12532-022-00223-3

Preprint: https://arxiv.org/abs/2006.05378

## BibTeX
```bibtex
@article{jalving2022graph,
  title   = {{A graph-based modeling abstraction for optimization: concepts and implementation in Plasmo.jl}},
  author  = {Jalving, Jordan and Shin, Sungho and Zavala, Victor M.},
  journal = {Mathematical Programming Computation},
  volume  = {14},
  number  = {4},
  pages   = {699--747},
  year    = {2022},
  doi     = {10.1007/s12532-022-00223-3}
}
```

## Authors
Jordan Jalving, Sungho Shin, and Victor M. Zavala
