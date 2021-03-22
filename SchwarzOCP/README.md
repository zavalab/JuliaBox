To reproduce the results, run the following script
```julia
julia --project=. -p 20 # make sure that your hardware has enough cores

] add https://github.com/sshin23/SimpleNLModels.jl
] add https://github.com/sshin23/SimpleSchwarz.jl
] add https://github.com/sshin23/SimpleADMM.jl
] instantiate

include("run_case_study.jl")
```
