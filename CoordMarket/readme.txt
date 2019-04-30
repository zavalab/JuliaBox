This folder provides supporting codes for the paper "Coordinated Management of Organic Waste and Derived Products".

* We have a Julia notebook file "Simple_Systems.ipynb" containing Julia code for illustrative case studies in the supporting information.

* We have a folder "Upper Yahara Case" for the main example in the paper. In this folder, we have four folders corresponding to four scenarios in the case study.

* In each case, we have three Julia scripts: "market_model.jl", "market_Run.jl", and "market_print.jl". One should run "market_Run.jl" first, this script will automatically read the "market_model.jl" script, establish the model, and solve the model. Then, the "market_print.jl" should be run in order to print out all the result files.

* If a sensitivity analysis on the VOS needs to be conducted (similar to the paper), one can change the lambda value in line 24 in "market_model.jl".

* We recommend use Julia 0.6.4 and Gurobi 8.1 to run all code files.


