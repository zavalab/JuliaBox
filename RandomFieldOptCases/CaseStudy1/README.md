# Dynamic Electricity Market
Here is the source-code for the simple dynamic optimal design case study that 
features random field uncertainty over time for market pricing. The entirety 
of this case study is self-contained in `dynamic_power.jl`. 

![power](power_network.png)

## Setting up the solver
This uses the Gurobi commercial solver which must be installed before setting up 
the Julia dependencies. Alternatively, one can modify the project file to use 
Ipopt instead (though this will take significantly longer to solve).

## Running it
To configure the required packages, we recommend creating a Julia environment 
using the included `Project.toml` file. Thus, we can configure the environment and 
run the case study via:
```julia
julia> cd("[INSERT_PATH_TO_FILES]/RandomFieldOptCases/CaseStudy1/")

julia> ]

(@v1.6) pkg> activate .

(CaseStudy1) pkg> instantiate

julia> include("dynamic_power.jl")
```
Note it will be slow the first time it is run as the packages are installed 
and precompiled. However, subsequent runs should be relatively quick.
