# Hierarchical Graphs

This repository contains scripts for building the 3-layer hierarchical optimization problem presented by [Atakan et al., 2022](https://doi-org.ezproxy.library.wisc.edu/10.1016/j.ejor.2021.12.042). We build this problem in Plasmo.jl and solve with Gurobi. These scripts were used in the manuscript "Hierarchical Graph Modeling for Multi-Scale Optimization of Power Systems". The repository is structured as follows

### datasets/
This directory contains the datasets from Atakan et al. (originally from the NREL-118 bus system from [Pena et al.](https://doi.org/10.1109/TPWRS.2017.2695963)). Some of these datasets have been processed through the files `HP_process_time_data.jl` and `data_interpolation.pynb`. 

### HP_low_reserves
This directory contains scripts for building the hierarchical optimization problem of Atakan et al. (with a few minor adjustments, as outlined in the manuscript above). We used two approaches for solving the hierarchical problem. First is to solve using a receding horizon approach where each subproblem is solved in series (32d_serial). The other approach was to solve an entire day as a monolithic problem (32d_monolith) This directory contains the following files or directories.

 * `layer_construction.jl` - contains functions for building the hierarchical structure. Also contains functions for linking subproblems between layers
 * `link_solutions.jl` - contains functions for linking *solutions* of previous problems to additional problems. Most of these functions are adaptations of functions in the `layer_construction.jl` script, where, instead of linking variables, we are linking the value of variables that have already been solved for. This is primarily used for the receding horizon approach, but the function `link_between_DA` links solutions from the previous days monolithic graph to the current day. 
 * `run_32d_serial.jl` and `run_32d_monolith.jl` - scripts for running all 32 days of the receding horizon or monolithic approaches. Data is saved in the `results/` directory. 
 * `test_up_down_times.jl` - contains code for testing whether the up and down time constraints were followed. These constraints proved the most difficult to model, so we wrote this script to test to ensure that the up and down time constraints were in deed being met in the results of the model. 
 * `32_day_data_analysis.jl` - contains code for extracting the solutions of the result files and putting in a dataframe with all 32 days of solutions
 * `get_objective_HA` - contains code for extracting the objective value of the optimization problem for both the receding horizon and monolithic problems. 
 * `results/` - directory containing subdirectories (`32d_serial` and `32d_monolith`) containing CSVs with solutions of the receding horizon and monolithic approches. Also contains 6 CSVs with data for the total number of generators committed, turned on, or turned off in the DAUC/STUC layers, and the amount of power overgenerated/curtailed, sent to the grid, and the amount of unmet load for each time point of the DAUC, STUC, and HAED layers. 
 * `gurobi_logs/` - directory containing the output logs of the monolithic problem

### HP_verylow_reserves
This directory contains almost identical scripts as the `HP_low_reserves` directory. The main difference is that the files for optimizing over the 32 days contain a different reserves requirement. In the work of Atakan et al., they give different reserve requirements for the optimization problem, modeled as a percentage of the demand. We used both the low reserves and very low reserves options given by Atakan et al. We only solved the very low reserves for the receding horizon problem, but it could be done for both approaches. 

### gephi_scripts
This directory contains scripts for extracting the graph structures as SimpleGraphs (from either LightGraphs.jl or Graphs.jl). The CSVs created from these scripts were used in the Gephi software to create visualizations of some of these problems. 

### Packages and Versions
To obtain the results in this work, we used the following packages and versions through Julia 1.7.3
 * Plasmo (v0.5.3)
 * Gurobi (v1.0.1)
 * DelimitedFiles
 * CSV (v0.10.9)
 * DataFrames (v1.5.0)
 * LightGraphs (v1.3.5)
 * Graphs (v1.8.0)
