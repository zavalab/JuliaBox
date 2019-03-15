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
- data/*.jld include data sets from pglib. https://github.com/power-grid-lib/pglib-opf
- The example.jl run the case study examples.

- The following MATLAB scripts generate the plots.
plotter.m

- To run the script, use the following command.
mpirun -np (number_of_processors) julia example.jl