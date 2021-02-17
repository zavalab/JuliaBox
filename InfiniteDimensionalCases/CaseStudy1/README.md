Here is the source-code for the stochastic optimal power flow case study that 
features different event constraints (e.g., joint chance constraints). The entirety 
of this case study is self-contained in `stochastic_optimal_powerflow.jl`. 

To configure the required packages, we recommend creating a Julia environment 
using the included `Project.toml` file. Thus, we can configure the environment and 
run the case study via:
```julia
julia> cd("[INSERT_PATH_TO_FILES]")

julia> ]

(@v1.5) pkg> activate .

julia> include("stochastic_optimal_powerflow.jl")
```
