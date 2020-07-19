This folder provides supporting codes for the paper "Multi-Site Supply Planning for Drug Products Under Uncertainty".
** We note that due to the proprietary nature of the data, the data files provided here are synthetic data files and are different the ones used for analysis in the paper. In order to facilitate comparison between the deterministic and stochastic solutions, the synthetic data are generated so that they produce trends similar to the inventory trends presented in the paper. 

* The folder "deterministic" contains the code and synthetic data files to solve the deterministic supply planning problem.

* The folder "stochastic" contains the code and synthetic data files to solve the stochastic supply planning problem.

* The folder "comparing_solutions_1000_scenarios" contains the code and synthetic data files to solve the stochastic supply planning problem.

* In each folder, we have six Julia scripts: "run_supply_planning_real.jl", "supply_model_backlog_st90.jl", "load_data.jl", "local_pre_plot.jl", "plot_results_new.jl", and "scenario_generator.jl". One should run "run_supply_planning_real.jl" first, this script will automatically read the "supply_model_backlog_st90.jl" and "load_data.jl" script, establish the model, and solve the model. Then, the "plot_results_new.jl" should be run in order to generate the inventory projection plots. The optional script of "scenario_generator.jl" is available if the interested reader wants to generate new scenarios for product yield, demand, and unplanned downtimes. 

* We recommend use Julia 0.6.4 and Gurobi 8.1 to run all code files.
