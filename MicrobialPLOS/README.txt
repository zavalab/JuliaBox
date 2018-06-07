This folder includes the Julia scripts used in manuscript "Scalable Nonlinear Programming Framework for Parameter Estimation in Dynamic Biological System Models" by Sungho Shin, Ophelia Venturelli, and Victor Zavala.

To run the julia scripts, the following Julia packages need to be installed:
JuMP, Plasmo, Ipopt, JLD.

The directions for the installation can be found in the following links:
https://github.com/JuliaOpt/JuMP.jl
https://github.com/jalving/Plasmo.jl
https://github.com/JuliaIO/JLD.jl
https://github.com/JuliaOpt/Ipopt.jl

The following is the explanation on each of the scripts.

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
plot_.m
plot_err_dist.m
plot_fitting.m
plot_hmap.m
plot_inference.m
run_all
