Here is the source-code for the pandemic optimal control problem where we
compare using CVaR against using an integral and/or a maximum measure for the 
objective function. The entirety of this case study is self-contained in 
`pandemic_control.jl`. 

To configure the required packages, we recommend creating a Julia environment 
using the included `Project.toml` file. Thus, we can configure the environment and 
run the case study via:
```julia
julia> cd("[INSERT_PATH_TO_FILES]")

julia> ]

(@v1.5) pkg> activate .

julia> include("pandemic_control.jl")
```
