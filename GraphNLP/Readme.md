# Gas Network Case Study Code

This repository contains the code for a stochastic gas network case study that we implemented in Julia. It is a two-stage program that maximized profit over different demand scenarios. The problem is modeled in Plasmo.jl and solved with MadNLP.jl. 

There are 3 Julia files

* modelfunctions_stoch.jl - This file contains functions to create subgraphs of components within the pipeline
* load_data_stoch.jl - This file contains a function to load in data from the JLD2 file "13_pipelines_200_scenarios.jld2" that contains data for each scenario. 
* run_gas_network.jl - This file calls the needed functions to create an OptiGraph containing subgraphs for each scenario. It also links first stage variables (compressor power) to a master node, and then it solves the optimization problem. 
