# Gas Network Case Study Code

This repository contains the code for a stochastic gas network case study that we implemented in Julia. It is a two-stage program that maximized profit over different demand scenarios. The problem is modeled in Plasmo.jl and solved with MadNLP.jl. 

There are 4 Julia files and 2 JLD2 files

* modelfunctions_stoch.jl - This file contains functions to create subgraphs of components within the pipeline
* load_data_stoch.jl - This file contains a function to load in data from the two JLD2 files that contain data for different scenarios
* run_gas_network.jl - This file calls the needed functions to create an OptiGraph containing subgraphs for each scenario. It also links first stage variables (compressor power) to a master node, and then it solves the optimization problem. 
* make_jld_stoch_150.jl - This file creates the JLD2 file with the data for each scenario.
* 13_pipelines_10_case_study.jld2 - This JLD2 file contains the data used to solve a 10 scenario case study
* 13_pipelines_150_case_study.jld2 - This JLD2 file contains the data used to scale the problem to up to 150 scenarios. 
