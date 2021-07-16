# Stochastic Optimal Pandemic Control
Here is the source-code for the pandemic optimal control problem where we
compare using CVaR against using an integral and/or a maximum measure for the 
objective function. The entirety of this case study is self-contained in 
`pandemic_control.jl`. 

![traj](covid_cvar_0_5.png)

## Running it
To configure the required packages, we recommend creating a Julia environment 
using the included `Project.toml` file. Thus, we can configure the environment and 
run the case study via:
```julia
julia> cd("[INSERT_PATH_TO_FILES]/InfiniteDimensionalCases/CaseStudy2/")

julia> ]

(@v1.6) pkg> activate .

(CaseStudy2) pkg> instantiate

julia> include("pandemic_control.jl")
```
Note it will be slow the first time it is run as the packages are installed 
and precompiled. However, subsequent runs should be relatively quick.
