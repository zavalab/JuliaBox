# Graph-Structured Nonlinear Optimization Code

This repository contains the code of two different nonlinear optimization problems for the manuscript "A Julia Framework for Graph-Structured Nonlinear Optimization." The problems are modeled in Plasmo.jl and solved in MadNLP.jl.

There are four folders containing code, and their contents are described below. 

## 1. stochastic_PID
This folder contains the code to model a stochastic PID controller problem in both JuMP.jl (StochPID_JuMP.jl) and Plasmo.jl (StochPID_Plasmo.jl). The code for Plasmo.jl also gives an example of how to manually repartition a graph-based model. Snippets of this code are contained in the manuscript mentioned above but are placed in this repository for convenience. The file StochPID_Plasmo.jl also contains code to plot the OptiGraphs partitioned in both scenario and time using the package PlasmoPlots.jl. The node locations generated from PlasmoPlots.jl was used to form some of the figures in the above manuscript. The code for forming these figures is given in the visualizations folder. There is also code in this file to create the adjacency matrices for the scenario and time partitioned OptiGraphs.

## 2. 10_scenario_gasnetwork
This folder contains code to solve a stochastic gas network case study. This is a two-stage stochastic program that maximized profit over different demand scenarios. The results of 10 scenarios are shown in the manuscript mentioned above. This folder contains code to run the model once and analyze the results. Note that the 10 scenario case study given in the manuscript uses the "13_pipelines_10_case_study.jld2" data. Details about the files in this folder are given below.

* modelfunctions_stoch.jl - This file contains functions to create subgraphs of components within the pipeline
* load_data_stoch.jl - This file contains a function to load in data from the two JLD2 files that contain data for different scenarios
* run_gas_network.jl - This file calls the needed functions to create an OptiGraph containing subgraphs for each scenario. It also links first stage variables (compressor power) to a master node, and then it solves the optimization problem. This file also contains code to plot a single scenario of the gasnetwork using PlasmoPlots.jl. The node locations generated from PlasmoPlots.jl was used to form some of the figures in the above manuscript. The code for forming these figures is given in the visualizations folder.
* make_jld_stoch_150.jl - This file creates the JLD2 file with the data for each scenario. The JLD2 files are presented here, but this file allows a user to more easily alter the jld2 data if they desire a different scenario.
* 13_pipelines_10_case_study.jld2 - This JLD2 file contains the data used to solve a 10 scenario case study
* 13_pipelines_150_case_study.jld2 - This JLD2 file contains the data used to scale the problem to up to 150 scenarios. 
* 10_scenario_results.csv - This file contains the optimal values for all the variables in the 10 scenarios of the problem presented in the manuscript above. The first column of the CSV gives the variable names and the 10 remaining columns are the optimal values for each scenario, where each column corresponds to a scenario
* 10_scenario_master_node.csv - This file contains the optimal compressor power (first stage variables) from the master node
* 10_scenario_analysis.jl - This file is an example script of how to work with the two CSV files of results. It contains scripts to get the linepack, amount of gas delivered, and compressor power from the above CSVs. It also contains the compressor and pipeline order for the data above. 

## 3. timed_gasnetwork
This folder contains the code used to time the gasnetwork problem at different problem scales. These scripts were used to obtain the solver times for MadNLP.jl operating with MA57, the Schur solver, and PardisoMKL. This folder also contains all of the output files for the data presented in the above manuscript. 

* modelfunctions_stoch.jl - This file contains functions to create subgraphs of components within the pipeline
* load_data_stoch.jl - This file contains a function to load in data from the two JLD2 files that contain data for different scenarios
* run_gas_network.jl - This file contains two funtions used to run the gas network model 4 times for a given solver and a given number of scenarios. The first of the 4 runs was omitted in our data analysis to avoid measuring compilation time. 
* 13_pipelines_150_case_study.jld2 - This JLD2 file contains the data used to scale the problem to up to 150 scenarios. 
* output_data.jl - This file contains the total solver time, linear solver time, and number of iterations for runs 2-4 from run_gas_network.jl tests. This same data is contained in the output files, but this file contains the data as arrays in Julia for convenience.
* output_files - This folder contains all of the runs for MadNLP.jl equipped with Ma57 running in both serial and parallel (30 threads), MadNLP.jl equipped with Schur solver running in parallel (30 threads), and MadNLP.jl equipped with PardisoMKL running in parallel (30 threads). Each of these runs was using a Plasmo.jl OptiGraph model of the gasnetwork problem. This folder also contains the runs for MadNLP.jl equipped with MA57 running in parallel (30 threads) for a JuMP.jl model rather than a Plasmo.jl OptiGraph. File names are "MadNLP_{run identifier}_{number of scenarios} {nx value} {run number}"

## 4. visualizations
This folder contains CSVs of node locations and a .jl file that was used to generate the figures for the above manuscript. The node locations were generated from PlasmoPlots (examples of the code for PlasmoPlots.jl can be seen in stochastic_PID/StochPID_Plasmo.jl and in 10_scenario_gasnetwork/run_gas_network.jl).

