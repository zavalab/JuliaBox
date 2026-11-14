# HierarchicalADMM

Hierarchical (ADMM-based) optimization architecture for large-scale power networks.

This repository includes the Julia/MATLAB scripts used in paper "A Hierarchical Optimization Architecture for Large-Scale Power Networks" by Sungho Shin (sungho.shin@wisc.edu), Philip Hart, Thomas Jahns, Victor M. Zavala* (victor.zavala@wisc.edu)

# Installation

To run the julia scripts, the following Julia packages need to be installed:

- JuMP.jl (v0.6.4) https://github.com/JuliaOpt/JuMP.jl
- Ipopt.jl (v0.3.0) https://github.com/JuliaOpt/Ipopt.jl
- JLD.jl (v0.8.3). https://github.com/JuliaIO/JLD.jl

Furthermore, to run the scripts with IPOPT with HSL linear solvers (e.g., MA57), one needs to install:

- IPOPT (v3.12.4) https://projects.coin-or.org/Ipopt
- HSL linear solvers http://www.hsl.rl.ac.uk/ipopt/

# Description
- multigrid.jl includes functions used in the example scripts.
- jlds/*.jld include data sets from pglib. https://github.com/power-grid-lib/pglib-opf
- The example.jl run the case study examples.

- The following MATLAB scripts generate the plots.
plotter.m

- To run the script, use the following command.
```shell
mpirun -np (number_of_processors) julia example.jl
```

## Reference
This code accompanies the paper:

Sungho Shin, Philip Hart, Thomas Jahns, and Victor M. Zavala. "A Hierarchical Optimization Architecture for Large-Scale Power Networks." *IEEE Transactions on Control of Network Systems* 6, no. 3 (2019): 1004-1014. https://doi.org/10.1109/tcns.2019.2906917

Preprint: https://arxiv.org/abs/2002.09796

## BibTeX
```bibtex
@article{shin2019hierarchical,
  title   = {{A Hierarchical Optimization Architecture for Large-Scale Power Networks}},
  author  = {Shin, Sungho and Hart, Philip and Jahns, Thomas and Zavala, Victor M.},
  journal = {IEEE Transactions on Control of Network Systems},
  volume  = {6},
  number  = {3},
  pages   = {1004--1014},
  year    = {2019},
  doi     = {10.1109/tcns.2019.2906917}
}
```

## Authors
Sungho Shin, Philip Hart, Thomas Jahns, and Victor M. Zavala
