This repository includes the Julia/MATLAB scripts used in manuscript "Scalable Nonlinear Programming Framework for Parameter Estimation in Dynamic Biological System Models" by Sungho Shin (sungho.shin@wisc.edu), Ophelia S. Venturelli (venturelli@wisc.edu), and Victor M. Zavala (victor.zavala@wisc.edu).

# Installation

To run the julia scripts, the following Julia packages need to be installed:

JuMP.jl (v0.6.4), Plasmo.jl (commit d00420598e406331f057a6eadd0ff7e54b691d71), Ipopt.jl (v0.3.0), JLD.jl (v0.8.3).

The directions for the installation can be found in the following links:

https://github.com/JuliaOpt/JuMP.jl
https://github.com/jalving/Plasmo.jl
https://github.com/JuliaIO/JLD.jl
https://github.com/JuliaOpt/Ipopt.jl

Furthermore, to run the scripts with PIPS-NLP and IPOPT with HSL linear solvers (e.g., MA57), one needs to install:

IPOPT (v3.12.4), PIPS-NLP (commit a1f347c157c742524cd90106c3c04c7abf2d0338).

The directions for the installation can be found in the following links:

https://projects.coin-or.org/Ipopt
https://github.com/Argonne-National-Laboratory/PIPS/tree/master/PIPS-NLP

# Description
- param.jl includes functions used in the example scripts.

- data.jld includes experimental/synthetic data sets.

- The following Julia scripts run the parameter estimation with different settings/data.
example_L1prior.jl
example_UQ.jl
example_cvar.jl
example_ipopt_P2.jl
example_mle.jl
example_pips_P2.jl
example_pips_S1.jl
example_pips_S2.jl
example_pips_S3.jl
example_pips_S4.jl
example_pips_bash
example_saturable.jl
example_standard.jl

- The following MATLAB scripts generate the plots.
plot_scalability.m
plot_err_dist.m
plot_fitting.m
plot_hmap.m
plot_inference.m

- run_all runs all the examples.
