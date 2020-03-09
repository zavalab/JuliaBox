This folder provides supporting codes for the paper "Valuing Economic Impact Reductions of Nutrient Pollution from Livestock Waste".

* The folder "sensitivity_analysis" contains the code and data files for different values of the economic impact/value of service (vos).

* The folder "GIS_data" contains the code and data files used to generate the maps of the Upper Yahara watershed region presented in the paper. 

* In each case, we have three Julia scripts: "market_model.jl", "market_Run.jl", and "market_print.jl". One should run "market_Run.jl" first, this script will automatically read the "market_model.jl" script, establish the model, and solve the model. Then, the "market_print.jl" should be run in order to print out all the result files.

* If a sensitivity analysis on the VOS needs to be conducted (similar to the paper), one can change the lambda value in line 26 in "market_model.jl".

* We recommend use Julia 0.6.4 and Gurobi 8.1 to run all code files for sensitivity analysis.
